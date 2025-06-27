#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <math.h>
#include <algorithm>
#include <unordered_set>
#include <bits/stdc++.h>

// 2D static array graph implementation
class Graph {
public:
    int n;             // number of vertices
    int m = 0;         // number of edges
    int blocks;        // number of 64-bit blocks per row
    std::vector<uint64_t> adjMatrix;
    std::vector<uint64_t> adjMatrixComplement;

    Graph(int numVertices)
        : n(numVertices),
          blocks((numVertices + 63) >> 6),
          adjMatrix(n * blocks, 0ULL),
          adjMatrixComplement(n * blocks, 0ULL)
    {
        uint64_t lastMask = ((blocks << 6) > n)
                          ? (~0ULL >> ((blocks << 6) - n))
                          : ~0ULL;
        for (int i = 0; i < n; ++i) {
            for (int b = 0; b < blocks; ++b)
                adjMatrixComplement[i*blocks + b] = ~0ULL;
            int bi = i >> 6, oi = i & 63;
            adjMatrixComplement[i*blocks + bi] &= ~(1ULL << oi);
            adjMatrixComplement[i*blocks + blocks - 1] &= lastMask;
        }
    }

    static Graph fromFile(const std::string &fileName) {
        std::ifstream file(fileName);
        if (!file) throw std::runtime_error("Cannot open " + fileName);
        int N; file >> N;
        Graph g(N);
        int u, v;
        while (file >> u >> v)
            g.addEdge(u, v);
        return g;
    }

    int numVertices() const { return n; }
    int numEdges() const    { return m; }

    inline bool hasEdge(int u, int v) const {
        return adjMatrix[u*blocks + (v>>6)] & (1ULL << (v&63));
    }
    inline bool hasCompEdge(int u, int v) const {
        return adjMatrixComplement[u*blocks + (v>>6)] & (1ULL << (v&63));
    }

    void addEdge(int u, int v) {
        if (u != v && u>=0 && v>=0 && u<n && v<n) {
            adjMatrix[u*blocks + (v>>6)] |= (1ULL << (v&63));
            adjMatrix[v*blocks + (u>>6)] |= (1ULL << (u&63));
            adjMatrixComplement[u*blocks + (v>>6)] &= ~(1ULL << (v&63));
            adjMatrixComplement[v*blocks + (u>>6)] &= ~(1ULL << (u&63));
            ++m;
        }
    }

    int degree(int u) const {
        int d = 0;
        for (int b = 0; b < blocks; ++b)
            d += __builtin_popcountll(adjMatrix[u*blocks + b]);
        return d;
    }

    const uint64_t* row(int u) const { return &adjMatrix[u*blocks]; }
    const uint64_t* complementRow(int u) const { return &adjMatrixComplement[u*blocks]; }
};

// ALGORITHM (START)

void bronKerbosch(const Graph &g,
                  const std::vector<int> &R,
                  const std::vector<int> &P,
                  const std::vector<int> &X,
                  const std::vector<bool> &active,
                  std::vector<int> &best)
{
    if (P.empty() && X.empty())
    {
        if (R.size() > best.size())
            best = R;
        return;
    }
    int pivot = -1, maxCnt = -1;
    for (int u : P)
    {
        if (!active[u])
            continue;
        const uint64_t *prow = g.complementRow(u);
        int cnt = 0;
        for (int v : P)
            if (active[v])
            {
                int b = v >> 6, o = v & 63;
                if (prow[b] & (1ULL << o))
                    ++cnt;
            }
        if (cnt > maxCnt)
        {
            maxCnt = cnt;
            pivot = u;
        }
    }
    const uint64_t *prow = (pivot >= 0 ? g.complementRow(pivot) : nullptr);
    std::vector<int> ext;
    for (int v : P)
    {
        if (!active[v])
            continue;
        if (pivot < 0 || !(prow[v >> 6] & (1ULL << (v & 63))))
            ext.push_back(v);
    }
    for (int v : ext)
    {
        std::vector<int> R2 = R, P2, X2;
        R2.push_back(v);
        const uint64_t *rowCv = g.complementRow(v);
        for (int w : P)
            if (active[w])
            {
                int b = w >> 6, o = w & 63;
                if (rowCv[b] & (1ULL << o))
                    P2.push_back(w);
            }
        for (int w : X)
            if (active[w])
            {
                int b = w >> 6, o = w & 63;
                if (rowCv[b] & (1ULL << o))
                    X2.push_back(w);
            }
        bronKerbosch(g, R2, P2, X2, active, best);
    }
}

static int greedyColourUB(const Graph &g,
                          const std::vector<int> &order,
                          const std::vector<bool> &active)
{
    int used = 0;
    std::vector<int> colour(g.n, -1);
    for (int v : order)
    {
        if (!active[v])
            continue;
        std::vector<char> forbidden(used, 0);
        const uint64_t *prow = g.complementRow(v);
        for (int u : order)
        {
            if (active[u] && colour[u] >= 0)
            {
                int b = u >> 6, o = u & 63;
                if (prow[b] & (1ULL << o))
                    forbidden[colour[u]] = 1;
            }
        }
        int c = 0;
        while (c < used && forbidden[c])
            ++c;
        if (c == used)
            ++used;
        colour[v] = c;
    }
    return used;
}

// include/exclude recursion
void branchAndBoundRec(const Graph &g,
                       const std::vector<int> &order,
                       std::vector<bool> &active,
                       std::vector<int> &curr,
                       std::vector<int> &best)
{
    int rem = 0;
    for (int v : order)
        if (active[v])
            ++rem;
    if ((int)curr.size() + rem <= (int)best.size())
        return;

    int ub = greedyColourUB(g, order, active);
    if ((int)curr.size() + ub <= (int)best.size())
        return;

    int v = -1;
    for (int u : order)
        if (active[u])
        {
            v = u;
            break;
        }
    if (v < 0)
    {
        if (curr.size() > best.size())
            best = curr;
        return;
    }

    std::vector<int> saved;
    saved.push_back(v);
    active[v] = false;
    const uint64_t *rowV = g.row(v);
    for (int u = 0; u < g.n; ++u)
        if (active[u])
        {
            int b = u >> 6, o = u & 63;
            if (rowV[b] & (1ULL << o))
            {
                active[u] = false;
                saved.push_back(u);
            }
        }
    curr.push_back(v);
    branchAndBoundRec(g, order, active, curr, best);
    curr.pop_back();
    for (int u : saved)
        active[u] = true;

    active[v] = false;
    branchAndBoundRec(g, order, active, curr, best);
    active[v] = true;
}
void reduceGraph(const Graph &g,
                 std::vector<bool> &active,
                 std::vector<int> &partial)
{
    bool changed = true;
    while (changed)
    {
        changed = false;
        for (int u = 0; u < g.n; ++u)
        {
            if (!active[u])
                continue;
            const uint64_t *row = g.row(u);
            std::vector<int> neigh;
            for (int v = 0; v < g.n; ++v)
                if (active[v])
                {
                    int b = v >> 6, o = v & 63;
                    if (row[b] & (1ULL << o))
                        neigh.push_back(v);
                }
            if (neigh.empty())
            {
                partial.push_back(u);
                active[u] = false;
                changed = true;
                break;
            }
            else if (neigh.size() == 1)
            {
                int v0 = neigh[0];
                partial.push_back(u);
                active[u] = false;
                active[v0] = false;
                changed = true;
                break;
            }
        }
    }
}

int branchAndBound(const Graph &g, std::vector<int> &independentSet)
{
    int n = g.n;
    std::vector<bool> active(n, true);
    std::vector<int> partial;
    reduceGraph(g, active, partial);

    std::vector<std::pair<int, int>> deg;
    for (int u = 0; u < n; ++u)
        if (active[u])
            deg.emplace_back(g.degree(u), u);
    std::sort(deg.begin(), deg.end(), [](auto &a, auto &b)
              { return a.first > b.first; });
    std::vector<int> order;
    for (auto &p : deg)
        order.push_back(p.second);

    std::vector<int> best = partial;
    for (int v : order)
    {
        bool ok = true;
        const uint64_t *rowV = g.row(v);
        for (int u : best)
        {
            int b = u >> 6, o = u & 63;
            if (rowV[b] & (1ULL << o))
            {
                ok = false;
                break;
            }
        }
        if (ok)
            best.push_back(v);
    }

    int rem = order.size(), edges = g.numEdges();
    double density = double(edges * 2) / (rem * (rem - 1));

    if (density >= 0.6)
    {
        std::vector<int> R = partial, P = order, X;
        bronKerbosch(g, R, P, X, active, best);
    }
    else
    {
        std::vector<int> current = partial;
        branchAndBoundRec(g, order, active, current, best);
    }

    independentSet = best;
    return best.size();
}

// ALGORITHM (END)

int main(int argc, char *argv[])
{
    std::string inputFile = argv[1];
    free(argv);

    std::unique_ptr<Graph> graph = std::make_unique<Graph>(Graph::fromFile(inputFile));

    std::vector<int> independentSet;

    auto start = std::chrono::high_resolution_clock::now();
    int size = branchAndBound(*graph, independentSet);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;

    std::sort(independentSet.begin(), independentSet.end());
    std::cout << "MIS size: " << size << "\n";
    std::cout << "Vertices: ";
    for (int v : independentSet)
        std::cout << v << " ";
    std::cout << "\n";
    std::cout << "Execution time: " << elapsed.count() << " seconds\n";

    return 0;
}