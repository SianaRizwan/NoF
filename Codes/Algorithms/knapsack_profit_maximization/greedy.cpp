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
#include <sys/resource.h>
#include <malloc.h>

using namespace std;
const int MAX_ITEMS = 10000;


void set_memory_limit(long mb_limit) {
    rlimit mem_limit{};
    mem_limit.rlim_cur = mb_limit * 1024 * 1024* 1024; // Convert MB to bytes
    mem_limit.rlim_max = mem_limit.rlim_cur;
    setrlimit(RLIMIT_AS, &mem_limit); // Set virtual memory limit
}


long getUsedMemoryKB() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss;
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







double greedyKnapsack(double weights[], double profits[], int n, double capacity, double &time_taken, long &mem_allocated)
{
    auto start = clock();

    int temp_sort[MAX_ITEMS];
    for (int i = 0; i < n; ++i)
        temp_sort[i] = i;

    // Sort profit-to-weight ratio (descending)
    sort(temp_sort, temp_sort + n, [&](int i, int j)
         { return (profits[i] / weights[i]) > (profits[j] / weights[j]); });

    // Greedy selection
    double total_weight = 0.0, total_profit = 0.0;
    int i = 0;
    while (i < n && total_weight < capacity)
    {
        int idx = temp_sort[i];
        if (weights[idx] <= (capacity - total_weight))
        {
            total_weight += weights[idx];
            total_profit += profits[idx];
        }
        ++i;
    }

    time_taken = 1000.0 * (clock() - start) / (double)CLOCKS_PER_SEC;
    mem_allocated = getUsedMemoryKB();

    return total_profit;
}






int main() {
    int ram = 256; //Set ram
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
        double profit = greedyKnapsack(weights, profits, n, instance.capacity, time_taken, peak_memory);
       
        
        outfile << n << "," << instance.capacity << "," << profit << "," << time_taken 
            << "," << ram << ",32," << peak_memory << "\n";
    }
    cout << "Finished processing " << instances.size() << " instances.\n";
    return 0;
}


