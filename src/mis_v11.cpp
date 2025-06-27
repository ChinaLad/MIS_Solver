#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <math.h>
#include <algorithm>
#include <unordered_set>
#include <bits/stdc++.h>

static inline int popcount64(uint64_t x) { return __builtin_popcountll(x); }

struct Bitset
{
    int blocks;           // number of 64-bit blocks
    const uint64_t *data; // pointer to start

    Bitset(int b = 0, const uint64_t *d = nullptr) : blocks(b), data(d) {}

    // count bits in active ∧ row(u)
    int intersectionCount(const Bitset &other) const
    {
        int cnt = 0;
        for (int b = 0; b < blocks; ++b)
            cnt += popcount64(data[b] & other.data[b]);
        return cnt;
    }
};

// 2D static array graph implementation
class Graph
{
private:
    int n;
    int m; // number of edges
    int blocks;

public:
    std::vector<uint64_t> adjMatrix;
    std::vector<uint64_t> adjMatrixComplement;
    int *degrees;

    Graph(int numVertices) : n(numVertices), blocks((numVertices + 63) >> 6), adjMatrix(n * blocks), adjMatrixComplement(n * blocks)
    {
        uint64_t allOnes = ~0ULL;
        for (int i = 0; i < n; ++i)
        {
            uint64_t *row = &adjMatrixComplement[i * blocks];
            for (int b = 0; b < blocks; ++b)
                adjMatrixComplement[i * blocks + b] = allOnes;
            // clear self-loop
            row[i >> 6] &= ~(1ULL << (i & 63));
            // clear extra bits
            int extra = (blocks << 6) - n;
            if (extra)
            {
                row[blocks - 1] &= (allOnes >> extra);
            }
        }
    }

    static Graph fromFile(const std::string &fileName)
    {
        std::ifstream file(fileName);
        if (!file)
            throw std::runtime_error("Cannot open " + fileName);
        int N;
        file >> N;
        Graph g(N);
        int u, v;
        while (file >> u >> v)
            g.addEdge(u, v);
        return g;
    }

    int numVertices() const { return n; }
    int numEdges() const { return m; }
    int numBlocks() const { return blocks; }

    bool hasEdge(int u, int v) const
    {
        return adjMatrix[u * blocks + (v >> 6)] & (1ULL << (v & 63));
    }

    void addEdge(int u, int v)
    {
        if (u != v && u >= 0 && v >= 0 && u < n && v < n)
        {
            adjMatrix[u * blocks + (v >> 6)] |= (1ULL << (v & 63));
            adjMatrix[v * blocks + (u >> 6)] |= (1ULL << (u & 63));
            adjMatrixComplement[u * blocks + (v >> 6)] &= ~(1ULL << (v & 63));
            adjMatrixComplement[v * blocks + (u >> 6)] &= ~(1ULL << (u & 63));
            degrees[u]++;
            degrees[v]++;
            m++;
        }
    }

    int degree(int u) const
    {
        return degrees[u];
    }

    // direct row-accessor for fast scans
    const uint64_t *row(int u) const { return &adjMatrix[u * blocks]; }

    Bitset rowBits(int u) const { return Bitset(blocks, (uint64_t *)&adjMatrix[u * blocks]); }

    const uint64_t *complementRow(int u) const { return &adjMatrixComplement[u * blocks]; }

    Bitset compRowBits(int u) const { return Bitset(blocks, (uint64_t *)&adjMatrixComplement[u * blocks]); }
};

// ALGORITHM (START)

void bronKerbosch(const Graph &g,
                  const std::vector<int> &R,
                  const std::vector<int> &P,
                  const std::vector<int> &X,
                  const std::vector<bool> &active,
                  std::vector<int> &best)
{
    // get complement matrix
    if (P.empty() && X.empty())
    {
        if (R.size() > best.size())
            best = R;
        return;
    }
    int pivot = -1, maxCnt = -1, n = g.numVertices(), blocks = g.numBlocks();
    // pivot select
    for (int u : P)
    {
        if (!active[u])
            continue;

        int cnt = Bitset(blocks, &g.adjMatrixComplement[u * blocks]).intersectionCount(Bitset(blocks, &g.adjMatrixComplement[0]));

        if (cnt > maxCnt)
        {
            maxCnt = cnt;
            pivot = u;
        }
    }

    Bitset pivotMask = (pivot > 0 ? Bitset(blocks, &g.adjMatrixComplement[pivot * blocks]) : Bitset());
    std::vector<int> ext;

    for (int v : P)
    {
        if (!active[v])
            continue;
        if (pivot < 0 || !(pivotMask.data[v >> 6] & (1ULL << (v & 63))))
            ext.push_back(v);
    }
    for (int v : ext)
    {
        std::vector<int> R2 = R;
        R2.push_back(v);
        std::vector<int> P2, X2;
        Bitset rowCv(blocks, &g.adjMatrixComplement[v * blocks]);
        for (int w : P)
            if (active[w] && (rowCv.data[w >> 6] & (1ULL << (w & 63))))
                P2.push_back(w);
        for (int w : X)
            if (active[w] && (rowCv.data[w >> 6] & (1ULL << (w & 63))))
                X2.push_back(w);
        bronKerbosch(g, R2, P2, X2, active, best);
    }
}

static int greedyColourUB(const Graph &g,
                          const std::vector<int> &order,
                          const std::vector<bool> &active)
{
    int n = g.numVertices(), blocks = g.numBlocks();
    std::vector<int> colour(n, -1);
    std::vector<std::vector<uint64_t>> colMask(64);
    int used = 0;

    // mask of active vertices as bitsets
    std::vector<uint64_t> activeMask(blocks);
    for (int v = 0; v < n; v++)
    {
        if (active[v])
        {
            activeMask[v >> 6] |= (1ULL << (v & 63));
        }
    }

    for (int v : order)
    {
        if (!active[v])
            continue;

        const uint64_t *rowC = g.complementRow(v);
        int c = 0;
        for (; c < used; c++)
        {
            bool conflict = false;
            for (int b = 0; b < blocks; b++)
            {
                if (rowC[b] & colMask[c][b])
                {
                    conflict = true;
                    break;
                }
            }
            if (!conflict)
                break;
        }
        if (c == used)
        {
            colMask.emplace_back(blocks);
            used++;
        }

        colMask[c][v >> 6] |= (1ULL << (v & 63));
    }
    return used;
}

// include/exclude recursion
void branchAndBoundRec(const Graph &g,
                       const std::vector<int> &order,
                       std::vector<bool> &active,
                       std::vector<int> &current,
                       std::vector<int> &best)
{
    int n = g.numVertices(), blocks = g.numBlocks();
    std::vector<uint64_t> activeMask(blocks);
    for (int v : order)
    {
        if (active[v])
            activeMask[v >> 6] |= (1ULL << (v & 63));
    }

    int rem = 0;
    for (int v : order)
        if (active[v])
            ++rem;

    // trivial bound
    if ((int)current.size() + rem <= (int)best.size())
        return;

    // colouring bound
    int ub = greedyColourUB(g, order, active);
    if ((int)current.size() + ub <= (int)best.size())
        return;

    // pick v
    int v = -1;
    for (int u : order)
        if (active[u])
        {
            v = u;
            break;
        }
    if (v < 0)
    {
        if (current.size() > best.size())
            best = current;
        return;
    }
    // include v
    std::vector<int> saved;
    saved.push_back(v);
    active[v] = false;
    const uint64_t *rowV = g.row(v);
    for (int b = 0; b < blocks; b++)
    {
        uint64_t bits = rowV[b] & activeMask[b];
        while (bits)
        {
            int u = __builtin_ctzll(bits);
            active[b * 64 + u] = false;
            saved.push_back(b * 64 + u);
            bits &= bits - 1;
        }
    }
    current.push_back(v);
    branchAndBoundRec(g, order, active, current, best);
    current.pop_back();
    for (int u : saved)
        active[u] = true;

    // exclude v
    active[v] = false;
    branchAndBoundRec(g, order, active, current, best);
    active[v] = true;
}

void reduceGraph(const Graph &g, std::vector<bool> &active, std::vector<int> &partial)
{
    int n = g.numVertices(), blocks = g.numBlocks();

    std::vector<uint64_t> activeMask(blocks);
    bool changed = true;
    while (changed)
    {
        changed = false;
        // update activeMask
        std::fill(activeMask.begin(), activeMask.end(), 0);
        for (int v = 0; v < n; v++)
            if (active[v])
                activeMask[v >> 6] |= 1ULL << (v & 63);

        for (int u = 0; u < n; ++u)
        {
            if (!active[u])
                continue;
            // degree in O
            uint64_t *row = (uint64_t *)&g.adjMatrix[u * blocks];
            int degU = 0;
            for (int b = 0; b < blocks; b++)
                degU += popcount64(row[b] & activeMask[b]);

            if (degU == 0)
            {
                partial.push_back(u);
                active[u] = false;
                changed = true;
                break;
            }
            else if (degU == 1)
            {
                // find the one neighbour
                for (int b = 0; b < blocks; b++)
                {
                    uint64_t bits = row[b] & activeMask[b];
                    if (bits)
                    {
                        int v = (b << 6) + __builtin_ctzll(bits);
                        partial.push_back(u);
                        active[u] = active[v] = false;
                        changed = true;
                        break;
                    }
                }
                if (changed)
                    break;
            }
        }
    }
}

int branchAndBound(const Graph &g, std::vector<int> &independentSet)
{
    int n = g.numVertices();
    std::vector<bool> active(n, true);
    std::vector<int> partial;
    reduceGraph(g, active, partial);

    std::vector<std::pair<int, int>> deg(n);
    for (int i = 0; i < n; ++i)
        deg[i] = {g.degree(i), i};

    // sort vertices by decreasing degree
    std::sort(deg.begin(), deg.end(), [](auto &a, auto &b)
              { return a.first > b.first; });
    std::vector<int> order(n);
    for (int i = 0; i < n; ++i)
        order[i] = deg[i].second;

    // greedy initial
    std::vector<int> best;
    for (int v : order)
    {
        bool ok = true;
        const uint64_t *row = g.row(v);
        for (int u : best)
            if (row[u >> 6] & (1ULL << (u & 63)))
            {
                ok = false;
                break;
            }
        if (ok)
            best.push_back(v);
    }

    int edges = g.numEdges();
    double density = (2.0 * edges) / (n * (n - 1));

    if (density >= 0.6)
    {
        std::vector<int> R = partial, P = order, X;
        bronKerbosch(g, R, P, X, active, best);
    }
    else
    {
        std::vector<int> curr = partial;
        branchAndBoundRec(g, order, active, curr, best);
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