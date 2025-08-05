#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <iomanip>
#include <stdexcept>
#include <random>
#include <sys/resource.h>
#include <malloc.h>
#include <limits>

using namespace std;

const int MAX_ITEMS = 10000;
const int MAX_POP = 100;       // max population size
const int MAX_GEN = 500;       // max generations
const int TOURNAMENT_K = 5;    // tournament size


void set_memory_limit(long mb_limit) {
    rlimit mem_limit{};
    mem_limit.rlim_cur = mb_limit * 1024 * 1024* 1024; // Convert GB to bytes
    mem_limit.rlim_max = mem_limit.rlim_cur;
    setrlimit(RLIMIT_AS, &mem_limit); // Set virtual memory limit
}

long getUsedMemoryKB() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss; // In kilobytes on Linux
}




struct Item {
    int id;
    double profit;
    double weight;
};

struct Instance {
    int num_elements;
    vector<Item> items;
    double capacity;
};

// Trims whitespace from a string
string trim(const string &s) {
    size_t start = s.find_first_not_of(" \t\n\r");
    size_t end = s.find_last_not_of(" \t\n\r");
    return (start == string::npos) ? "" : s.substr(start, end - start + 1);
}

vector<double> parseList(const string& str) {
    vector<double> result;
    string cleaned = str;

    cleaned.erase(remove(cleaned.begin(), cleaned.end(), '['), cleaned.end());
    cleaned.erase(remove(cleaned.begin(), cleaned.end(), ']'), cleaned.end());
    cleaned.erase(remove(cleaned.begin(), cleaned.end(), '"'), cleaned.end());  

    stringstream ss(cleaned);
    string item;
    while (getline(ss, item, ',')) {
        item = trim(item);
        if (!item.empty()) {
            result.push_back(stod(item));
        }
    }
    return result;
}


vector<Instance> parseCSV(const string& file_path) {
    vector<Instance> instances;
    ifstream file(file_path);
    string line;

    if (!file.is_open()) throw runtime_error("File open failed");

    getline(file, line); // Skip header
    int line_num = 1;

    while (getline(file, line)) {
        ++line_num;

        vector<string> row;
        string token;
        bool in_brackets = false;
        stringstream field;
        for (char c : line) {
            if (c == '[') in_brackets = true;
            if (c == ']') in_brackets = false;

            if (c == ',' && !in_brackets) {
                row.push_back(field.str());
                field.str("");
                field.clear();
            } else {
                field << c;
            }
        }
        row.push_back(field.str());

        if (row.size() != 4) {
            cerr << "Skipping line " << line_num << ": invalid column count\n";
            continue;
        }

        try {
            int num_elements = stoi(trim(row[0]));
            vector<double> weights = parseList(row[1]);
            vector<double> prices = parseList(row[2]);
            double capacity = stod(trim(row[3]));

            if ((int)weights.size() != num_elements || (int)prices.size() != num_elements) {
                cerr << "Skipping line " << line_num << ": mismatch in element count\n";
                continue;
            }

            vector<Item> items;
            for (int i = 0; i < num_elements; ++i) {
                items.push_back({i + 1, prices[i], weights[i]});
            }

            instances.push_back({num_elements, items, capacity});
        } catch (const exception& e) {
            cerr << "Skipping line " << line_num << ": parse error: " << e.what() << "\n";
        }
    }

    return instances;
}






double fitness(int individual[], const double profits[], const double weights[], int n, double capacity) {
    double total_profit = 0.0, total_weight = 0.0;
    for (int i = 0; i < n; ++i) {
        if (individual[i]) {
            total_profit += profits[i];
            total_weight += weights[i];
        }
    }

    // Penalize if weight < capacity (infeasible)
    if (total_weight + 1e-6 < capacity)
        return total_profit + 1e6 * (capacity - total_weight); // Large penalty

    return total_profit;  // Lower is better
}

// Create random bitstring 
void createIndividual(int indiv[], int n, mt19937& rng) {
  uniform_int_distribution<> bit(0,1);
  for(int i=0;i<n;++i)
    indiv[i] = bit(rng);
}

// Mutation by flipping bits
void mutate(int indiv[], int n, double mutation_rate, mt19937& rng) {
  uniform_real_distribution<> dist(0.0,1.0);
  for(int i=0;i<n;++i)
    if(dist(rng) < mutation_rate)
      indiv[i] = 1 - indiv[i];
}

// One-point crossover
void crossover(const int p1[], const int p2[], int child[], int n, mt19937& rng) {
  int point = uniform_int_distribution<>(1, n-1)(rng);
  for(int i=0;i<n;++i)
    child[i] = (i<point ? p1[i] : p2[i]);
  // child may be under- or over-weight; fitness will penalize under-weight
}


int tournamentSelection(int population[][MAX_ITEMS], const double profits[], const double weights[], int n, int pop_size, double capacity, mt19937& rng) {
    double best_fit = numeric_limits<double>::infinity();
    int best_idx = -1;
    for (int i = 0; i < TOURNAMENT_K; ++i) {
        int idx = uniform_int_distribution<>(0, pop_size - 1)(rng);
        double fit = fitness(population[idx], profits, weights, n, capacity);
        if (fit < best_fit) {
            best_fit = fit;
            best_idx = idx;
        }
    }
    return best_idx;
}

double geneticKnapsack(const double weights[], const double profits[], int n, double capacity, double& time_taken, long& mem_allocated, int population_size = 100,
    int generations = 300, double mutation_rate = 0.05, double timeout = 300.0) {

    clock_t start = clock();
    mt19937 rng(42);

    
    static int population[MAX_POP][MAX_ITEMS];
    static int next_gen [MAX_POP][MAX_ITEMS];

    // Initialize population with random bit‐strings
    for (int i = 0; i < population_size; ++i) {
        createIndividual(population[i], n, rng);
    }

    // Evolve
    for (int gen = 0; gen < generations; ++gen) {
        time_taken = (clock() - start) / (double)CLOCKS_PER_SEC;
        if (time_taken > timeout) {
            mem_allocated = getUsedMemoryKB();
            return numeric_limits<double>::infinity();
        }

        // Evaluate fitness of each individual
        static double fitness_scores[MAX_POP];
        static int    indices       [MAX_POP];
        for (int i = 0; i < population_size; ++i) {
            fitness_scores[i] = fitness(population[i], profits, weights, n, capacity);
            indices[i]        = i;
        }

        // Sort indices by increasing fitness ( minimize profit)
        sort(indices, indices + population_size,
             [&](int a, int b) {
                 return fitness_scores[a] < fitness_scores[b];
             });

        // Elitism: carry forward best half
        for (int i = 0; i < population_size / 2; ++i) {
            int src = indices[i];
            for (int j = 0; j < n; ++j)
                next_gen[i][j] = population[src][j];
        }

        // Fill rest by selection, crossover, mutation
        for (int i = population_size / 2; i < population_size; ++i) {
            int p1 = tournamentSelection(
                population, profits, weights, n,
                population_size, capacity, rng
            );
            int p2 = tournamentSelection(
                population, profits, weights, n,
                population_size, capacity, rng
            );

            crossover(population[p1], population[p2],
                      next_gen[i], n, rng);
            mutate(next_gen[i], n, mutation_rate, rng);
        }

        // Swap in next generation
        for (int i = 0; i < population_size; ++i)
            for (int j = 0; j < n; ++j)
                population[i][j] = next_gen[i][j];
    }

    //  Return best (min) profit among feasible solutions
    double best = numeric_limits<double>::infinity();
    for (int i = 0; i < population_size; ++i) {
        double f = fitness(population[i], profits, weights, n, capacity);
        if (f < best) best = f;
    }

    time_taken    = (clock() - start) / (double)CLOCKS_PER_SEC;
    mem_allocated = getUsedMemoryKB();
    return best;
}





int main() {
    int ram = 256;
    set_memory_limit(ram); 
    string input_file = "dataset_200.csv"; //Specify the input dataset filepath
    string output_file = "output.csv"; //Specify the output filepath
    auto instances = parseCSV(input_file);

    ofstream outfile(output_file);
    outfile << fixed << setprecision(6);
    outfile << "number_of_elements,capacity,total_profit,solution_time,ram,cpu_cores,peak_memory\n";

    for (const auto& instance : instances) {
        int n = instance.num_elements;
        double weights[MAX_ITEMS], profits[MAX_ITEMS];

        for (int i = 0; i < n; ++i) {
            weights[i] = instance.items[i].weight;
            profits[i] = instance.items[i].profit;
        }

        double time_taken;
        long peak_memory = 0;
        double ga_profit = geneticKnapsack(weights, profits, n, instance.capacity, time_taken, peak_memory);

      

        outfile << instance.num_elements << "," << instance.capacity << "," << ga_profit << "," << time_taken 
            << "," << ram << ",32," << peak_memory << "\n";
    }
    cout << "Finished processing " << instances.size() << " instances.\n";
    return 0;
}