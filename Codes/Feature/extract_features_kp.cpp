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

using namespace std;

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

    if (!file.is_open()) throw runtime_error("Cannot open file");

    getline(file, line); // Skip header
    int line_num = 1;

    while (getline(file, line)) {
        ++line_num;

        vector<string> row;
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
            vector<double> profits = parseList(row[2]);
            double capacity = stod(trim(row[3]));

            if ((int)weights.size() != num_elements || (int)profits.size() != num_elements) {
                cerr << "Skipping line " << line_num << ": mismatch in element count\n";
                continue;
            }

            vector<Item> items;
            for (int i = 0; i < num_elements; ++i) {
                items.push_back({i + 1, profits[i], weights[i]});
            }

            instances.push_back({num_elements, items, capacity});
        } catch (const exception& e) {
            cerr << "Skipping line " << line_num << ": parse error: " << e.what() << "\n";
        }
    }

    return instances;
}

// helpers
 double mean(const vector<double>& vec) {
    return vec.empty() ? 0.0 : accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
}

double median(vector<double> vec) {
    if (vec.empty()) return 0.0;
    sort(vec.begin(), vec.end());
    size_t n = vec.size();
    return (n % 2 == 0) ? (vec[n/2 - 1] + vec[n/2]) / 2.0 : vec[n/2];
}

double stdev(const vector<double>& vec, double m) {
    if (vec.size() < 2) return 0.0;
    double sum = 0.0;
    for (double v : vec) sum += (v - m) * (v - m);
    return sqrt(sum / (vec.size() - 1));
}

// Feature extraction
map<string, double> computeStatistics(const Instance& instance) {
    vector<double> weights, profits;
    for (const auto& item : instance.items) {
        weights.push_back(item.weight);
        profits.push_back(item.profit);
    }

    double max_w = *max_element(weights.begin(), weights.end());
    double min_w = *min_element(weights.begin(), weights.end());
    double mean_w = mean(weights);
    double median_w = median(weights);
    double std_w = stdev(weights, mean_w);

    double max_p = *max_element(profits.begin(), profits.end());
    double min_p = *min_element(profits.begin(), profits.end());
    double mean_p = mean(profits);
    double median_p = median(profits);
    double std_p = stdev(profits, mean_p);

    // Pearson correlation
    double correlation = 0.0;
    if (weights.size() > 1) {
        double cov = 0.0;
        for (size_t i = 0; i < weights.size(); ++i)
            cov += (weights[i] - mean_w) * (profits[i] - mean_p);
        cov /= weights.size();
        correlation = (std_w != 0.0 && std_p != 0.0) ? cov / (std_w * std_p) : 0.0;
    }

    double total_profit = accumulate(profits.begin(), profits.end(), 0.0);
    return {
        {"max_weight", max_w}, {"min_weight", min_w}, {"mean_weight", mean_w},
        {"median_weight", median_w}, {"std_weight", std_w}, {"weight_range", max_w - min_w},
        {"max_profit", max_p}, {"min_profit", min_p},
        {"mean_profit", mean_p}, {"median_profit", median_p}, {"std_profit", std_p}, {"profit_range", max_p - min_p},
        {"renting_ratio", total_profit != 0.0 ? instance.capacity / total_profit : 0.0},
        {"mean_weight_profit_ratio", mean_p != 0.0 ? mean_w / mean_p : 0.0},
        {"median_weight_profit_ratio", median_p != 0.0? median_w / median_p : 0.0},
        {"capacity_mean_weight_ratio", mean_w != 0.0 ? instance.capacity / mean_w : 0.0},
        {"capacity_median_weight_ratio", median_w != 0.0 ? instance.capacity / median_w : 0.0},
        {"capacity_std_weight_ratio", std_w != 0.0 ? instance.capacity / std_w : 0.0},
        {"std_weight_profit_ratio", std_p != 0.0 ? std_w / std_p : 0.0},
        {"weight_profit_correlation", correlation},
    };
}

int main() {
    string input_file = "dataset_200.csv"; //Specify path to dataset
    string output_file = "./results_cpp/features_kp.csv"; //Specify path to output file
    auto instances = parseCSV(input_file);

    ofstream outfile(output_file);
    outfile << fixed << setprecision(6);
    outfile << "number_of_elements,capacity,max_weight,min_weight,mean_weight,median_weight,std_weight,weight_range,"
            << "max_profit,min_profit,mean_profit,median_profit,std_profit,profit_range,renting_ratio,"
            << "mean_weight_profit_ratio,median_weight_profit_ratio,capacity_mean_weight_ratio,"
            << "capacity_median_weight_ratio,capacity_std_weight_ratio,std_weight_profit_ratio,"
            << "weight_profit_correlation\n";

    for (const auto& instance : instances) {
        auto stats = computeStatistics(instance);
        outfile << instance.num_elements << "," << instance.capacity;
        outfile << "," << stats["max_weight"] << "," << stats["min_weight"];
        outfile << "," << stats["mean_weight"] << "," << stats["median_weight"];
        outfile << "," << stats["std_weight"] << "," << stats["weight_range"];
        outfile << "," << stats["max_profit"] << "," << stats["min_profit"];
        outfile << "," << stats["mean_profit"] << "," << stats["median_profit"];
        outfile << "," << stats["std_profit"] << "," << stats["profit_range"];
        outfile << "," << stats["renting_ratio"];
        outfile << "," << stats["mean_weight_profit_ratio"];
        outfile << "," << stats["median_weight_profit_ratio"];
        outfile << "," << stats["capacity_mean_weight_ratio"];
        outfile << "," << stats["capacity_median_weight_ratio"];
        outfile << "," << stats["capacity_std_weight_ratio"];
        outfile << "," << stats["std_weight_profit_ratio"];
        outfile << "," << stats["weight_profit_correlation"];
        outfile << "\n";
    }

    cout << "Finished processing " << instances.size() << " instances.\n";
    cout << "Features written to " << output_file << "\n";
    return 0;
}
