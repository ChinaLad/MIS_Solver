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
    int m; // number of edges
    bool **adjMatrix;
    bool **adjMatrixComplement;
    int *degrees;

public:
    Graph(int numVertices) : n(numVertices)
    {
        adjMatrix = new bool *[n];
        adjMatrixComplement = new bool *[n];

        degrees = new int[n];
        std::fill(degrees, degrees + n, 0);

        for (int i = 0; i < n; ++i)
        {
            adjMatrix[i] = new bool[n];
            adjMatrixComplement[i] = new bool[n];
            std::fill(adjMatrix[i], adjMatrix[i] + n, false);
            std::fill(adjMatrixComplement[i], adjMatrixComplement[i] + n, true);
            adjMatrixComplement[i][i] = false; // no self-loops in complement
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

    bool hasEdge(int u, int v) const
    {
        return adjMatrix[u][v];
    }

    void addEdge(int u, int v)
    {
        if (u != v && u >= 0 && v >= 0 && u < n && v < n)
        {
            adjMatrix[u][v] = adjMatrix[v][u] = true;
            adjMatrixComplement[u][v] = adjMatrixComplement[v][u] = false;
            degrees[u]++;
            degrees[v]++;
            m++;
        }
    }

    void deleteEdge(int u, int v)
    {
        if (u != v && u >= 0 && v >= 0 && u < n && v < n)
        {
            adjMatrix[u][v] = adjMatrix[v][u] = false;
            adjMatrixComplement[u][v] = adjMatrixComplement[v][u] = true;
            degrees[u]--;
            degrees[v]--;
            m--;
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

    const bool *compRow(int u) const { return adjMatrixComplement[u]; }

    bool **complement() const
    {
        return adjMatrixComplement;
    }
};

// ALGORITHM (START)
/**
 * Algo v6 + reduction to maximum clique problem with complement graph as input.
 * Uses Bron-Kerbosch algorithm with pivoting.
 */

void bronKerbosch(const Graph &g,
                  std::vector<int> R,
                  std::vector<int> P,
                  std::vector<int> X,
                  const std::vector<bool> &active,
                  std::vector<int> &best)
{
    if (P.empty() && X.empty())
    {
        if (R.size() > best.size())
        {
            std::cout << "‣ new clique size " << R.size()
                      << "   R = { ";
            for (int v : R)
                std::cout << v << ' ';
            std::cout << "}\n";
            best = R;
        }
        return;
    }
    // choose pivot u from P∪X maximizing |P ∩ N_comp(u)|
    int pivot = -1, maxCnt = -1;
    auto consider = [&](int u)
    {
        int cnt = 0;
        const bool *comp = g.compRow(u);
        for (int v : P)
            if (active[v] && comp[v])
                ++cnt;
        if (cnt > maxCnt)
        {
            maxCnt = cnt;
            pivot = u;
        }
    };
    for (int u : P)
        if (active[u])
            consider(u);
    for (int u : X)
        if (active[u])
            consider(u);

    // build ext = P \ N_comp(pivot)
    std::vector<int> ext;
    if (pivot < 0)
    {
        ext = P;
    }
    else
    {
        const bool *pComp = g.compRow(pivot);
        for (int v : P)
            if (active[v] && (pivot < 0 || !pComp[v]))
                ext.push_back(v);
    }

    std::cout << "pivot=" << pivot
          << "  |P|=" << P.size()
          << "  |X|=" << X.size()
          << "  ext=" << ext.size() << '\n';
    // pivoting recursion
    for (int v : ext)
    {
        R.push_back(v);
        // P2 = P ∩ N_comp(v), X2 = X ∩ N_comp(v)
        std::vector<int> P2, X2;
        const bool *cComp = g.compRow(v);
        for (int w : P)
            if (active[w] && cComp[w])
                P2.push_back(w);
        for (int w : X)
            if (active[w] && cComp[w])
                X2.push_back(w);
        bronKerbosch(g, R, P2, X2, active, best);
        // move v from P to X
        P.erase(std::find(P.begin(), P.end(), v));
        X.push_back(v);
        R.pop_back();
    }
}

// include/exclude recursion
void branchAndBoundRec(const Graph &g,
                       const std::vector<int> &order,
                       std::vector<bool> &active,
                       std::vector<int> &current,
                       std::vector<int> &best)
{
    int n = g.numVertices();
    int rem = 0;
    for (int v : order)
        if (active[v])
            ++rem;
    if ((int)current.size() + rem <= (int)best.size())
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
        if (current.size() > best.size())
            best = current;
        return;
    }
    // include v
    std::vector<int> saved;
    active[v] = false;
    saved.push_back(v);
    const bool *rowV = g.row(v);
    for (int u = 0; u < n; ++u)
        if (active[u] && rowV[u])
        {
            active[u] = false;
            saved.push_back(u);
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

int branchAndBound(const Graph &g, std::vector<int> &independentSet)
{
    int n = g.numVertices();
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

    int edges = g.numEdges();
    double density = (2.0 * edges) / (n * (n - 1));

    std::vector<bool> active(n, true);
    std::vector<int> current;
    if (density >= 0.5)
    {
        std::vector<int> R, P = order, X;
        std::cout << "Using Bron-Kerbosch with pivoting\n";
        bronKerbosch(g, R, P, X, active, best);
    }
    else
    {
        std::cout << "Using Branch and Bound\n";
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