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
#include <queue>
#include <sys/resource.h>
#include <malloc.h>

using namespace std;
const size_t MAX_QUEUE_SIZE = 1000000;
const double BOUND_EPSILON = 1e-3;  // Ignore nodes with marginal bound improvement

void set_memory_limit(long mb_limit) {
    rlimit mem_limit{};
    mem_limit.rlim_cur = mb_limit * 1024 * 1024 * 1024;  // Convert GB to bytes
    mem_limit.rlim_max = mem_limit.rlim_cur;
    setrlimit(RLIMIT_AS, &mem_limit);  // Set virtual memory limit
}

long getUsedMemoryKB() {
    rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss;  // In kilobytes on Linux
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
    return (start == string::npos) ? string() : s.substr(start, end - start + 1);
}

vector<double> parseList(const string &str) {
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

    if (!file.is_open()) throw runtime_error("Cannot open file");

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


// Branch-and-Bound node
struct BBNode {
    double upper_bound;
    double profit;
    double weight;
    int level;
    bool operator<(const BBNode &other) const {
        return upper_bound < other.upper_bound;  // max-heap (by upper bound)
    }
};

// Function to solve knapsack with Branch-and-Bound
double solveKnapsackBB(const vector<Item>& raw_items, double capacity, double& time_taken, long &mem_allocated, double timeout = 300.0) {
    auto start_time = clock();
    if (capacity <= 0 || raw_items.empty()) {
        time_taken = (clock() - start_time) / (double)CLOCKS_PER_SEC;
        mem_allocated = getUsedMemoryKB();
        return 0.0;
    }

    int n = raw_items.size();
    int cap = static_cast<int>(capacity);

    auto items = raw_items;
    sort(items.begin(), items.end(), [](const Item& a, const Item& b) {
        return (a.profit / a.weight) > (b.profit / b.weight);
    });

    auto calculate_upper_bound = [&](int level, double profit, double weight) {
        if (weight > capacity) return 0.0;
        double bound = profit;
        double total_weight = weight;
        int i = level + 1;
        while (i < n && total_weight + items[i].weight <= capacity) {
            total_weight += items[i].weight;
            bound += items[i].profit;
            ++i;
        }
        if (i < n) {
            bound += (capacity - total_weight) * (items[i].profit / items[i].weight);
        }
        return bound;
    };

    priority_queue<BBNode> pq;
    pq.push({INFINITY, 0.0, 0.0, -1});

    double best_profit = 0.0;

    while (!pq.empty()) {
        if ((clock() - start_time) / (double)CLOCKS_PER_SEC > timeout) break;

        BBNode node = pq.top();
        pq.pop();

        if (node.upper_bound <= best_profit + BOUND_EPSILON) continue;

        int level = node.level + 1;
        if (level >= n) continue;

        const Item& current = items[level];

        // Include current item
        double new_profit = node.profit + current.profit;
        double new_weight = node.weight + current.weight;

        if (new_weight <= capacity && new_profit > best_profit) {
            best_profit = new_profit;
        }

        double new_bound = calculate_upper_bound(level, new_profit, new_weight);
        if (new_bound > best_profit + BOUND_EPSILON && pq.size() < MAX_QUEUE_SIZE) {
            pq.push({new_bound, new_profit, new_weight, level});
        }

        // Exclude current item
        new_bound = calculate_upper_bound(level, node.profit, node.weight);
        if (new_bound > best_profit + BOUND_EPSILON && pq.size() < MAX_QUEUE_SIZE) {
            pq.push({new_bound, node.profit, node.weight, level});
        }
    }

    time_taken = (clock() - start_time) / (double)CLOCKS_PER_SEC;
    mem_allocated = getUsedMemoryKB();
    return best_profit;
}



int main() {
    int ram = 256; //Specify RAM
    set_memory_limit(ram);
    string input_file = "dataset_200.csv"; //Specify the path of the input dataset
    string output_file = "./output.csv"; //Specify the path of the output file
    auto instances = parseCSV(input_file);

    ofstream outfile(output_file);
    outfile << fixed << setprecision(6);
    outfile << "number_of_elements,capacity,profit,solution_time,ram,cpu_cores,peak_memory\n";

    for (const auto &instance : instances) {
        double time_taken;
        long peak_memory = 0;
        double profit = solveKnapsackBB(instance.items, instance.capacity, time_taken, peak_memory);

        outfile << instance.num_elements << "," << instance.capacity << ","
                << profit << "," << time_taken << "," << ram << ",32," << peak_memory << "\n";
    }

    cout << "Processed" << instances.size() << " instances.\n";
    return 0;
}
