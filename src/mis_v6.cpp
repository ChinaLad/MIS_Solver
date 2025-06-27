#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <math.h>
#include <algorithm>
#include <unordered_set>
#include <bits/stdc++.h>

// 2D static array graph implementation
class Graph
{
private:
    int n;
    bool **adjMatrix;
    int *degrees;

public:
    Graph(int numVertices) : n(numVertices)
    {
        adjMatrix = new bool *[n];
        degrees = new int[n];
        std::fill(degrees, degrees + n, 0);

        for (int i = 0; i < n; ++i)
        {
            adjMatrix[i] = new bool[n];
            std::fill(adjMatrix[i], adjMatrix[i] + n, false);
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

    bool hasEdge(int u, int v) const
    {
        return adjMatrix[u][v];
    }

    void addEdge(int u, int v)
    {
        if (u != v && u >= 0 && v >= 0 && u < n && v < n)
        {
            adjMatrix[u][v] = adjMatrix[v][u] = true;
            degrees[u]++;
            degrees[v]++;
        }
    }

    void deleteEdge(int u, int v)
    {
        if (u != v && u >= 0 && v >= 0 && u < n && v < n)
        {
            adjMatrix[u][v] = adjMatrix[v][u] = false;
            degrees[u]--;
            degrees[v]--;
        }
    }

    int degree(int u) const
    {
        return degrees[u];
    }

    std::vector<int> getNeighbours(int u) const
    {
        std::vector<int> out;
        for (int v = 0; v < n; ++v)
            if (adjMatrix[u][v])
                out.push_back(v);
        return out;
    }

    // direct row-accessor for fast scans
    const bool *row(int u) const { return adjMatrix[u]; }
};

bool isValidIndependentSet(const Graph &g, const std::vector<int> &set)
{
    for (size_t i = 0; i < set.size(); ++i)
    {
        for (size_t j = i + 1; j < set.size(); ++j)
        {
            if (g.hasEdge(set[i], set[j]))
            {
                return false; // Found an edge between two vertices in the set
            }
        }
    }
    return true; // No edges found, valid independent set
}

// ALGORITHM (START)
/**
 *   Algo v5 + better branching that excludes neighbors of included vertices.
 */

void branchAndBoundRec(const Graph &g, const std::vector<int> &order, std::vector<bool> &active, std::vector<int> &current, std::vector<int> &best)
{
    int n = order.size();
    // count active vertices
    int rem = 0;
    for (int v : order)
        if (active[v])
            ++rem;
    // bound: even if include all remaining, cannot beat best
    if ((int)current.size() + rem <= (int)best.size())
        return;
    // find next active v
    int v = -1;
    for (int u : order)
        if (active[u])
        {
            v = u;
            break;
        }
    if (v == -1)
    {
        // no more candidates
        if (current.size() > best.size())
            best = current;
        return;
    }
    // include v and exclude its neighbors
    {
        // record neighbors to deactivate
        std::vector<int> changed;
        // deactivate v
        active[v] = false;
        changed.push_back(v);
        // deactivate all neighbors of v
        const bool *rowV = g.row(v);
        for (int u = 0; u < n; ++u)
        {
            if (active[u] && rowV[u])
            {
                active[u] = false;
                changed.push_back(u);
            }
        }
        current.push_back(v);
        branchAndBoundRec(g, order, active, current, best);
        current.pop_back();
        // restore
        for (int u : changed)
            active[u] = true;
    }

    // exclude v
    active[v] = false;
    branchAndBoundRec(g, order, active, current, best);
    active[v] = true;
}

int branchAndBound(const Graph &g, std::vector<int> &independentSet)
{
    int n = g.numVertices();
    // compute degrees
    std::vector<std::pair<int, int>> deg(n);
    for (int i = 0; i < n; ++i)
        deg[i] = {g.degree(i), i};

    // sort vertices by decreasing degree
    std::sort(deg.begin(), deg.end(), [](auto &a, auto &b)
              { return a.first > b.first; });
    std::vector<int> order(n);
    for (int i = 0; i < n; ++i)
        order[i] = deg[i].second;

    // greedy initial solution to seed best lower bound
    std::vector<int> best;
    for (int v : order)
    {
        bool ok = true;
        const bool *rowV = g.row(v);
        for (int u : best)
        {
            if (rowV[u])
            {
                ok = false;
                break;
            }
        }
        if (ok)
        {
            best.push_back(v);
        }
    }
    
    // recursive branch & bound
    std::vector<bool> active(n, true);
    std::vector<int> current;
    branchAndBoundRec(g, order, active, current, best);

    independentSet = std::move(best);
    return independentSet.size();
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