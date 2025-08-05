#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <string>
#include <limits>
#include <algorithm>
#include <sys/resource.h>
#include <filesystem>

using namespace std;
using namespace std::filesystem;


struct City {
    int id;
    double x, y;
};

void set_memory_limit(long gb_limit) {
    rlimit mem_limit{};
    mem_limit.rlim_cur = gb_limit * 1024LL * 1024LL * 1024LL;
    mem_limit.rlim_max = mem_limit.rlim_cur;
    setrlimit(RLIMIT_AS, &mem_limit);
}

long getUsedMemoryKB() {
    rusage usage{};
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss;
}

// n×n distance matrix
using DistMatrix = vector<vector<double>>;
DistMatrix buildDistMatrix(const vector<City>& cities) {
    int n = cities.size();
    DistMatrix D(n, vector<double>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            double d = hypot(cities[i].x - cities[j].x,
                             cities[i].y - cities[j].y);
            D[i][j] = D[j][i] = d;
        }
    }
    return D;
}

// Nearest‑Neighbor
double nearestNeighbor(const DistMatrix& D, vector<int>& tour) {
    int n = D.size();
    vector<bool> used(n, false);
    tour.clear();
    tour.push_back(0);
    used[0] = true;
    double cost = 0;
    for (int step = 1; step < n; ++step) {
        int cur = tour.back(), nxt = -1;
        double best = numeric_limits<double>::infinity();
        for (int j = 0; j < n; ++j) {
            if (!used[j] && D[cur][j] < best) {
                best = D[cur][j];
                nxt = j;
            }
        }
        tour.push_back(nxt);
        used[nxt] = true;
        cost += best;
    }
    cost += D[tour.back()][0];
    return cost;
}

int main() {
    int ramGB    = 256;
    int cpuCores = 32;

    set_memory_limit(ramGB);

    path inputDir = "tsp";
    path outPath   = path("tsp_cpp") / "nearest_neighbor" /
                         ("NN_" + to_string(cpuCores) +
                          "_"  + to_string(ramGB) + ".csv");
    create_directories(outPath.parent_path());

    ofstream ofs(outPath);
    ofs << "tsp_file,ram_gb,cpu_cores,total_cost,solution_time_ms,peak_memory_kb\n";
    ofs << fixed << setprecision(6);

    int count = 0;
    for (auto &ent : directory_iterator(inputDir)) {
        if (ent.path().extension() != ".tsp") continue;

        // read cities
        vector<City> cities;
        {
            ifstream in(ent.path());
            string line;
            // skip until coordinate section
            while (getline(in, line) &&
                   line.find("NODE_COORD_SECTION") == string::npos)
                ;
            // read until EOF
            while (getline(in, line) &&
                   line.find("EOF") == string::npos) {
                istringstream iss(line);
                City c;
                if (iss >> c.id >> c.x >> c.y)
                    cities.push_back(c);
            }
        }

        if (cities.size() < 2) {
            cerr << "Skipping \"" << ent.path().filename()
                      << "\": only " << cities.size() << " city read.\n";
            continue;
        }

        
        DistMatrix D = buildDistMatrix(cities);
        vector<int> tour;
        double bestCost;

        clock_t start = clock();
        bestCost = nearestNeighbor(D, tour);
        clock_t end   = clock();
        double timeMs = 1000.0 * (end - start) / CLOCKS_PER_SEC;
        if (timeMs > 300000.0) bestCost = 0.0;

        long peakKB = getUsedMemoryKB();

        ofs << ent.path().filename().string() << ","
            << ramGB << "," << cpuCores << ","
            << bestCost << "," << timeMs << "," << peakKB << "\n";
        ofs.flush();
        ++count;
    }

    cout << "Finished processing " << count << " instances.\n";
    return 0;
}
