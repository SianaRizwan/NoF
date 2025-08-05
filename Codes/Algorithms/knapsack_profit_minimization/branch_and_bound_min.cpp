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
#include <limits>

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
    double lower_bound;
    double profit;
    double weight;
    int level;
    bool operator<(const BBNode& other) const {
        return lower_bound > other.lower_bound;  // min‐heap
    }
};



double greedy_integer(const vector<Item>& items, double capacity) {
    double w = 0.0, p = 0.0;
    for (const auto& it : items) {
        w += it.weight;
        p += it.profit;
        if (w >= capacity) {
            // We've reached or exceeded capacity with whole items
            return p;
        }
    }
    // No integer fill reached capacity
    return numeric_limits<double>::infinity();
}

// Greedy fractional solver to get an initial finite best_profit
double greedy_upper_bound(const vector<Item>& items, double capacity) {
    double w = 0.0, p = 0.0;
    for (auto& it : items) {
        if (w + it.weight >= capacity) {
            double take = capacity - w;
            p += take * (it.profit / it.weight);
            return p;
        }
        w += it.weight;
        p += it.profit;
    }
    return numeric_limits<double>::infinity();
}


double solveKnapsackBB(const vector<Item>& raw_items,
                       double capacity,
                       double& time_taken,
                       long& mem_allocated,
                       double timeout = 300.0) {
    clock_t start_time = clock();

    if (capacity <= 0 || raw_items.empty()) {
        time_taken = (clock() - start_time) / (double)CLOCKS_PER_SEC;
        mem_allocated = getUsedMemoryKB();
        return 0.0;
    }

    int n = raw_items.size();
    auto items = raw_items;

    // Sort by increasing profit/weight
    sort(items.begin(), items.end(),
              [](const Item& a, const Item& b) {
                  return (a.profit / a.weight) < (b.profit / b.weight);
              });

    // Build prefix sums for weights & profits
    vector<double> cumWeight(n+1, 0.0), cumProfit(n+1, 0.0);
    for (int i = 0; i < n; ++i) {
        cumWeight[i+1] = cumWeight[i] + items[i].weight;
        cumProfit[i+1] = cumProfit[i] + items[i].profit;
    }

    // Initialize best_profit via greedy integer fill
    double best_profit = greedy_upper_bound(items, capacity);

    // Fast bound: binary search on prefix sums
    auto calculate_lower_bound = [&](int level, double profit, double weight) {
        if (weight >= capacity) return profit;
        double rem = capacity - weight;
        int lo = level + 1, hi = n;
        while (lo < hi) {
            int mid = (lo + hi + 1) / 2;
            if (cumWeight[mid] - cumWeight[level+1] <= rem) lo = mid;
            else hi = mid - 1;
        }
        double bound = profit + (cumProfit[lo] - cumProfit[level+1]);
        if (lo < n && (cumWeight[lo] - cumWeight[level+1]) < rem) {
            double fraction = (rem - (cumWeight[lo] - cumWeight[level+1]))
                              / items[lo].weight;
            bound += fraction * items[lo].profit;
        }
        return bound;
    };

    //Push root node
    priority_queue<BBNode> pq;
    double root_bound = calculate_lower_bound(-1, 0.0, 0.0);
    pq.push({ root_bound, 0.0, 0.0, -1 });

    // Branch‐and‐bound loop
    while (!pq.empty()) {
        double elapsed = (clock() - start_time) / (double)CLOCKS_PER_SEC;
        if (elapsed > timeout) break;

        BBNode node = pq.top(); pq.pop();
        if (node.lower_bound >= best_profit - 1e-9) continue;

        int level = node.level + 1;
        if (level >= n) continue;

        const Item& cur = items[level];

        // include
        double p_inc = node.profit + cur.profit;
        double w_inc = node.weight + cur.weight;
        if (w_inc >= capacity && p_inc < best_profit) {
            best_profit = p_inc;
        }
        double b_inc = calculate_lower_bound(level, p_inc, w_inc);
        if (b_inc < best_profit - 1e-9 && pq.size() < MAX_QUEUE_SIZE) {
            pq.push({ b_inc, p_inc, w_inc, level });
        }

        // exclude
        double b_exc = calculate_lower_bound(level, node.profit, node.weight);
        if (b_exc < best_profit - 1e-9 && pq.size() < MAX_QUEUE_SIZE) {
            pq.push({ b_exc, node.profit, node.weight, level });
        }
    }

    time_taken    = (clock() - start_time) / (double)CLOCKS_PER_SEC;
    mem_allocated = getUsedMemoryKB();
    return (best_profit == numeric_limits<double>::infinity())
           ? 0.0
           : best_profit;
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
