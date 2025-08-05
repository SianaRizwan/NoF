#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <map>
#include <set>
#include <filesystem>
#include <queue>
#include <limits>
#include <utility>

using namespace std;
using namespace fs;
namespace fs = filesystem;

typedef pair<double, double> Point;

struct Node { int id; double x, y; };

struct Stats {
    double min = 0;
    double max = 0;
    double mean = 0;
    double median = 0;
    double stddev = 0;
    double cv = 0;
};

// Summary statistics
Stats computeStats(vector<double>& v) {
    Stats s;
    if (v.empty()) return s;
    sort(v.begin(), v.end());
    s.min = v.front();
    s.max = v.back();
    s.mean = accumulate(v.begin(), v.end(), 0.0) / v.size();
    if (v.size() == 1) {
        s.median = s.mean;
        s.stddev = 0;
        s.cv = 0;
        return s;
    }
    if (v.size() % 2 == 0)
        s.median = (v[v.size()/2 - 1] + v[v.size()/2]) / 2.0;
    else
        s.median = v[v.size()/2];
    double sumsq = 0;
    for (double d : v) sumsq += (d - s.mean) * (d - s.mean);
    s.stddev = sqrt(sumsq / (v.size() - 1));
    s.cv = (s.mean != 0 ? s.stddev / s.mean : 0);
    return s;
}

//Parse .tsp file
bool parseTSP(const string& path, vector<Node>& nodes) {
    ifstream in(path);
    if (!in.is_open()) return false;
    string line;
    bool inSection = false;
    while (getline(in, line)) {
        if (line.rfind("NODE_COORD_SECTION", 0) == 0) { inSection = true; continue; }
        if (!inSection) continue;
        if (line == "EOF" || line.empty()) break;
        istringstream iss(line);
        Node node;
        if (!(iss >> node.id >> node.x >> node.y)) continue;
        nodes.push_back(node);
    }
    return !nodes.empty();
}

// Compute all pairwise distances
void computeDistances(const vector<Node>& nodes,
                      vector<vector<double>>& D,
                      vector<double>& distList) {
    int N = nodes.size();
    D.assign(N, vector<double>(N, 0));
    distList.clear();
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double d = hypot(nodes[i].x - nodes[j].x,
                                  nodes[i].y - nodes[j].y);
            D[i][j] = D[j][i] = d;
            distList.push_back(d);
        }
    }
}

// Prim's MST
void primMST(const vector<vector<double>>& D,
             vector<int>& parent,
             vector<double>& edgeWeights) {
    int N = D.size();
    parent.assign(N, -1);
    edgeWeights.assign(N, numeric_limits<double>::infinity());
    vector<bool> used(N, false);
    edgeWeights[0] = 0;
    for (int it = 0; it < N; ++it) {
        int u = -1;
        double best = numeric_limits<double>::infinity();
        for (int v = 0; v < N; ++v) {
            if (!used[v] && edgeWeights[v] < best) {
                best = edgeWeights[v];
                u = v;
            }
        }
        if (u < 0) break;
        used[u] = true;
        for (int v = 0; v < N; ++v) {
            if (!used[v] && D[u][v] < edgeWeights[v]) {
                parent[v] = u;
                edgeWeights[v] = D[u][v];
            }
        }
    }
    // remove the dummy root edge
    vector<double> temp;
    for (int i = 1; i < N; ++i) temp.push_back(edgeWeights[i]);
    edgeWeights.swap(temp);
}

// Compute depths via BFS
vector<int> computeDepths(int root, const vector<int>& parent) {
    int N = parent.size();
    vector<vector<int>> adj(N);
    for (int i = 0; i < N; ++i) {
        if (parent[i] >= 0) {
            adj[i].push_back(parent[i]);
            adj[parent[i]].push_back(i);
        }
    }
    vector<int> depth(N, 0);
    vector<bool> vis(N, false);
    queue<int> q;
    vis[root] = true;
    q.push(root);
    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int v : adj[u]) {
            if (!vis[v]) {
                vis[v] = true;
                depth[v] = depth[u] + 1;
                q.push(v);
            }
        }
    }
    return depth;
}

int main() {
    
    string folder = "tsp_200_data";
    path folderPath(folder);
    
    path outPath = folderPath / "tsp_features.csv";
    ofstream out(outPath.string());
  
    out << "instance_filename,"
           "edge_distance_min,edge_distance_max,edge_distance_mean,edge_distance_median,"
           "edge_distance_stddev,edge_distance_prop_below_mean,edge_distance_frac_distinct,"
           "expected_random_tour_length,edge_cost_mode_count,edge_cost_mode_peak_freq,"
           "edge_cost_mode_mean,"
           "cluster_count_eps0.01,cluster_mean_dist_centroid_eps0.01,"
           "cluster_count_eps0.05,cluster_mean_dist_centroid_eps0.05,"
           "cluster_count_eps0.1,cluster_mean_dist_centroid_eps0.1,"
           "nnd_min_rel,nnd_max_rel,nnd_mean_rel,nnd_median_rel,nnd_stddev_rel,nnd_cv_rel,"
           "centroid_x,centroid_y,centroid_dist_min,centroid_dist_mean,centroid_dist_max,"
           "mst_depth_min,mst_depth_max,mst_depth_mean,mst_depth_median,mst_depth_stddev,"
           "mst_edge_min,mst_edge_max,mst_edge_mean,mst_edge_median,mst_edge_stddev,"
           "mst_total_length_norm,"
           "angle_min,angle_max,angle_mean,angle_median,angle_stddev,"
           "convex_hull_area,convex_hull_frac_nodes"
        << endl;

    for (auto& p : directory_iterator(folderPath)) {
        if (p.path().extension() != ".tsp") continue;
        vector<Node> nodes;
        if (!parseTSP(p.path().string(), nodes)) continue;
        int N = nodes.size();

        // distances
        vector<vector<double>> D;
        vector<double> distList;
        computeDistances(nodes, D, distList);
        Stats ds = computeStats(distList);
        int cntBelow = count_if(distList.begin(), distList.end(), [&](double d){ return d < ds.mean; });
        double propBelow = distList.empty() ? 0 : (double)cntBelow / distList.size();
        set<double> distinct(distList.begin(), distList.end());
        double fracDistinct = distinct.empty() ? 0 : (double)distinct.size() / distList.size();
        double sumDist = accumulate(distList.begin(), distList.end(), 0.0);
        double expTour = (N > 1) ? sumDist * 2.0 / (N - 1) : 0;
        // mode statistics
        map<double,int> freq;
        for (double d : distList) freq[round(d*1e6)/1e6]++;
        int peakFreq = 0;
        for (auto& kv : freq) peakFreq = max(peakFreq, kv.second);
        vector<double> modes;
        for (auto& kv : freq) if (kv.second == peakFreq) modes.push_back(kv.first);
        int modeCount = modes.size();
        double modeMean = modeCount>0 ? accumulate(modes.begin(), modes.end(), 0.0)/modeCount : 0;

        // clustering via DBSCAN (minPts=3) 
        auto runDBSCAN = [&](double epsFrac, int minPts) {
            double epsDist = epsFrac * ds.max;
            vector<int> labels(N, -1);
            int C = 0;
            for (int i = 0; i < N; ++i) {
                if (labels[i] != -1) continue;
                // neighborhood
                vector<int> neighbors;
                for (int j = 0; j < N; ++j)
                    if (D[i][j] <= epsDist) neighbors.push_back(j);
                if (neighbors.size() < (size_t)minPts) {
                    labels[i] = -2;
                    continue;
                }
                // new cluster
                labels[i] = C;
                queue<int> q;
                for (int nb: neighbors) if (nb != i) q.push(nb);
                while (!q.empty()) {
                    int cur = q.front(); q.pop();
                    if (labels[cur] == -2) labels[cur] = C;
                    if (labels[cur] != -1) continue;
                    labels[cur] = C;
                    // expand
                    vector<int> curNeigh;
                    for (int k = 0; k < N; ++k)
                        if (D[cur][k] <= epsDist) curNeigh.push_back(k);
                    if (curNeigh.size() >= (size_t)minPts) {
                        for (int nb2: curNeigh)
                            if (labels[nb2] == -1) q.push(nb2);
                    }
                }
                C++;
            }
            if (C == 0) return make_pair(0, 0.0);
            // compute centroids
            vector<Point> cents(C,{0,0});
            vector<int> cnts(C,0);
            for (int i = 0; i < N; ++i) {
                if (labels[i] >= 0) {
                    cents[labels[i]].first += nodes[i].x;
                    cents[labels[i]].second += nodes[i].y;
                    cnts[labels[i]]++;
                }
            }
            for (int c = 0; c < C; ++c) {
                cents[c].first /= cnts[c];
                cents[c].second /= cnts[c];
            }
            // mean distance to centroid
            double sumd = 0; int tot = 0;
            for (int i = 0; i < N; ++i) {
                if (labels[i] >= 0) {
                    auto [cx2, cy2] = cents[labels[i]];
                    sumd += hypot(nodes[i].x - cx2, nodes[i].y - cy2);
                    tot++;
                }
            }
            double meanDist = tot>0 ? sumd/tot : 0.0;
            return make_pair(C, meanDist);
        };
        auto c1 = runDBSCAN(0.01, 3);
        auto c2 = runDBSCAN(0.05, 3);
        auto c3 = runDBSCAN(0.1, 3);

        // nearest neighbor distances
        vector<double> nnd;
        for (int i = 0; i < N; ++i) {
            double m = numeric_limits<double>::infinity();
            for (int j = 0; j < N; ++j) if (i != j) m = min(m, D[i][j]);
            nnd.push_back(m / (ds.max>0 ? ds.max : 1));
        }
        Stats nnds = computeStats(nnd);

        // centroid distances
        double cx = 0, cy = 0;
        for (auto& nd : nodes) { cx += nd.x; cy += nd.y; }
        cx /= N; cy /= N;
        vector<double> cDist;
        for (auto& nd: nodes) cDist.push_back(hypot(nd.x - cx, nd.y - cy));
        Stats cds = computeStats(cDist);

        // MST
        vector<int> parent;
        vector<double> mEdges;
        primMST(D, parent, mEdges);
        Stats mes = computeStats(mEdges);
        double normM = (sumDist>0) ? accumulate(mEdges.begin(), mEdges.end(), 0.0) / sumDist : 0;

        // MST depths
        auto depths = computeDepths(0, parent);
        vector<double> depths_d(depths.begin(), depths.end());
        Stats ds_depth = computeStats(depths_d);

        // angles
        vector<double> angles;
        for (int i = 0; i < N; ++i) {
            vector<pair<double,int>> nbr;
            for (int j = 0; j < N; ++j) if (i != j) nbr.emplace_back(D[i][j], j);
            sort(nbr.begin(), nbr.end());
            int a = nbr[0].second, b = nbr[1].second;
            double ux = nodes[a].x - nodes[i].x, uy = nodes[a].y - nodes[i].y;
            double vx = nodes[b].x - nodes[i].x, vy = nodes[b].y - nodes[i].y;
            double dot = ux*vx + uy*vy;
            double denom = hypot(ux,uy) * hypot(vx,vy);
            double cosang = (denom > 0) ? dot/denom : 1.0;
            cosang = min(1.0, max(-1.0, cosang));
            angles.push_back(acos(cosang));
        }
        Stats angs = computeStats(angles);

        // convex hull
        vector<Point> pts;
        for (auto& nd: nodes) pts.emplace_back(nd.x, nd.y);
        // compute convex hull (Monotone chain)
        auto cross = [](const Point& a, const Point& b, const Point& c){
            return (b.first - a.first)*(c.second - a.second) - (b.second - a.second)*(c.first - a.first);
        };
        sort(pts.begin(), pts.end());
        vector<Point> H;
        for (int pass = 0; pass < 2; ++pass) {
            int start = H.size();
            for (auto& pt: pts) {
                while ((int)H.size() >= start+2 && cross(H[H.size()-2], H.back(), pt) <= 0)
                    H.pop_back();
                H.push_back(pt);
            }
            H.pop_back();
            reverse(pts.begin(), pts.end());
        }
        double area = 0;
        for (size_t i = 1; i+1 < H.size(); ++i)
            area += (H[i].first - H[0].first)*(H[i+1].second - H[0].second)
                  - (H[i].second - H[0].second)*(H[i+1].first - H[0].first);
        area = abs(area)/2.0;
        double fracHull = H.empty() ? 0 : (double)H.size() / N;

     
        out << p.path().filename().string() << ","
            // edge stats
            << ds.min << "," << ds.max << "," << ds.mean << "," << ds.median << "," << ds.stddev << ","
            << propBelow << "," << fracDistinct << "," << expTour << ","
            // mode
            << modeCount << "," << peakFreq << "," << modeMean << ","
            // clusters
            << c1.first << "," << c1.second << ","
            << c2.first << "," << c2.second << ","
            << c3.first << "," << c3.second << ","
            // nnd
            << nnds.min << "," << nnds.max << "," << nnds.mean << "," << nnds.median << "," << nnds.stddev << "," << nnds.cv << ","
            // centroid
            << cx << "," << cy << "," << cds.min << "," << cds.mean << "," << cds.max << ","
            // MST depth
            << ds_depth.min << "," << ds_depth.max << "," << ds_depth.mean << "," << ds_depth.median << "," << ds_depth.stddev << ","
            // MST edge
            << mes.min << "," << mes.max << "," << mes.mean << "," << mes.median << "," << mes.stddev << ","
            << normM << ","
            // angles
            << angs.min << "," << angs.max << "," << angs.mean << "," << angs.median << "," << angs.stddev << ","
            // convex hull
            << area << "," << fracHull
            << endl;
    }
    out.close();
    cout << "Feature extraction completed" << endl;
    return 0;
}