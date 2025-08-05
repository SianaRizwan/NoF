#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <map>
#include <numeric>
#include <iomanip>
#include <stdexcept>
#include <sys/resource.h>
#include <malloc.h>


using namespace std;


const int MAX_ITEMS = 10000;
const int MAX_CAPACITY = 1000000; // Maximum capacity for the knapsack
double SCALE = 1.0;


void set_memory_limit(long mb_limit) {
    rlimit mem_limit{};
    mem_limit.rlim_cur = mb_limit * 1024 * 1024* 1024; // Convert MB to bytes
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







double dpKnapsack(double weights[], double profits[], int n, int  scaled_cap, double &time_taken, long  &mem_allocated)
{
    clock_t start = clock();
    const double INF = numeric_limits<double>::infinity();

    // dp[w] = min profit to reach exactly w,
    // and dp[scaled_cap] also aggregates all w >= scaled_cap
    vector<double> dp(scaled_cap+1, INF);
    dp[0] = 0.0;

    for(int i=0;i<n;++i){
        int w = static_cast<int>(weights[i]*SCALE + 0.5);
        double p = profits[i];
        for(int j=scaled_cap;j>=0;--j){
            if(dp[j]==INF) continue;
            int j2 = min(j+w, scaled_cap);
            dp[j2] = min(dp[j2], dp[j] + p);
        }
        time_taken = (clock()-start)/(double)CLOCKS_PER_SEC;
        if(time_taken>300.0){
            mem_allocated = getUsedMemoryKB();
            return INF;
        }
    }

    mem_allocated = getUsedMemoryKB();
    return dp[scaled_cap];
}






int main() {
    int ram = 256;
    set_memory_limit(ram); 
    string input_file = "dataset_200.csv"; //Specify the input dataset filepath
    string output_file = "output.csv"; //Specify the output dataset filepath
    auto instances = parseCSV(input_file);
    

    ofstream outfile(output_file);
    outfile << fixed << setprecision(6);
    outfile << "number_of_elements,capacity,total_profit,solution_time,ram,cpu_cores,peak_memory\n";

   for (const auto& instance : instances) {
    int n = instance.num_elements;
    double weights[MAX_ITEMS], profits[MAX_ITEMS];
    double max_w = 0.0;   
   

    for (int i = 0; i < n; ++i) {
        weights[i] = instance.items[i].weight;
        profits[i] = instance.items[i].profit;
        max_w = max(max_w, weights[i]);
    }

    double instance_max = max(instance.capacity, max_w);
    double local_scale = instance_max>0 ? double(MAX_CAPACITY)/instance_max: 1.0;
    SCALE = local_scale;
    int scaled_cap = int(instance.capacity * SCALE + 0.5);
    double time_taken;
    long peak_memory = 0;
    
       
    double dp_profit = dpKnapsack(weights, profits, n, scaled_cap, time_taken, peak_memory);
    outfile << instance.num_elements << "," << instance.capacity << "," << dp_profit << "," << time_taken 
            << "," << ram << ",32," << peak_memory << "\n";
}
    cout << "Finished processing " << instances.size() << " instances.\n";
    return 0;
}