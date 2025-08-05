#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <string>
#include <limits>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <sys/resource.h>
#include <filesystem>

using namespace std;
using namespace std::filesystem;

// City record
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

// Nearest‑Neighbor initializer
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

// 2‑opt move encoding for tabu
static long long encodeEdge(int a, int b) {
    if (a > b) swap(a, b);
    return ((long long)a << 32) | b;
}

// Tabu‑Search with randomized neighbor sampling
vector<int> tabuSearch(const DistMatrix& D,
                       double& bestCost,
                       int maxIter = 1000,
                       int tabuTenure = 50,
                       int neighSize = 100)
{
    int n = D.size();
    mt19937_64 rng(random_device{}());
    uniform_int_distribution<int> pick(1, n-1);

    vector<int> currTour;
    bestCost = nearestNeighbor(D, currTour);
    vector<int> bestTour = currTour;
    double currCost = bestCost;

    unordered_map<long long,int> tabu;

    for (int iter = 0; iter < maxIter; ++iter) {
        double bestDelta = numeric_limits<double>::infinity();
        int bi = -1, bj = -1;
        for (int k = 0; k < neighSize; ++k) {
            int i = pick(rng), j = pick(rng);
            if (i >= j) continue;
            int u = currTour[i-1], v = currTour[i];
            int x = currTour[j],   y = currTour[(j+1) % n];
            double delta = D[u][x] + D[v][y] - D[u][v] - D[x][y];
            long long moveCode = encodeEdge(v, x);
            bool isTabu = tabu.count(moveCode) && tabu[moveCode] > iter;
            if ((delta < bestDelta) && (!isTabu || currCost + delta < bestCost)) {
                bestDelta = delta; bi = i; bj = j;
            }
        }
        if (bi < 0) break;
        reverse(currTour.begin()+bi, currTour.begin()+bj+1);
        currCost += bestDelta;
        tabu[ encodeEdge(currTour[bi], currTour[bj]) ] = iter + tabuTenure;
        if (currCost < bestCost) {
            bestCost = currCost;
            bestTour = currTour;
        }
    }
    return bestTour;
}

int main() {
  
    int ramGB    = 256; 
    int cpuCores = 32;   

    set_memory_limit(ramGB);

    path inputDir = "tsp";
    path outPath   = path("tsp_cpp") / "tabu_search" /
                         ("TS_" + to_string(cpuCores) +
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
        double bestCost;

        
        clock_t start = clock();
        auto tour = tabuSearch(D, bestCost);
        clock_t end   = clock();
        double timeMs = 1000.0 * (end - start) / CLOCKS_PER_SEC;
        if (timeMs > 300000.0) bestCost = 0.0;

        long peakKB = getUsedMemoryKB();

        ofs << ent.path().filename().string() << ","
            << ramGB << "," << cpuCores << ","
            << bestCost << "," << timeMs << "," << peakKB << "\n";
        ofs.flush();
    }

    cout << "Finished "<< distance(directory_iterator(inputDir), directory_iterator{}) << " instances.\n";
    return 0;
}
