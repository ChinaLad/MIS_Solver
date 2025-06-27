#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <math.h>
#include <algorithm>
#include <unordered_set>
#include <bits/stdc++.h>

// 2D dynamic array graph implementation
class Graph
{
private:
    int n;
    std::vector<std::vector<bool>> adjMatrix;

public:
    Graph(int numVertices) : n(numVertices), adjMatrix(numVertices, std::vector<bool>(numVertices, false)) {}

    static Graph fromFile(const std::string &fileName)
    {
        std::ifstream file(fileName);
        if (!file.is_open())
        {
            std::cerr << "Error opening file: " << fileName << std::endl;
        }

        int numVertices;
        file >> numVertices;
        Graph g(numVertices);

        int u, v;
        while (file >> u >> v)
        {
            g.addEdge(u, v);
        }

        return g;
    }

    int numVertices() const
    {
        return n;
    }

    bool hasEdge(int u, int v) const
    {
        return adjMatrix[u][v];
    }

    void addEdge(int u, int v)
    {
        if (0 <= u && u < n && 0 <= v && v < n && u != v)
        {
            adjMatrix[u][v] = true;
            adjMatrix[v][u] = true;
        }
    }

    void deleteEdge(int u, int v)
    {
        if (0 <= u && u < n && 0 <= v && v < n && u != v)
        {
            adjMatrix[u][v] = false;
            adjMatrix[v][u] = false;
        }
    }

    std::vector<int> getNeighbours(int u) const
    {
        std::vector<int> neighbours;
        for (int v = 0; v < n; v++)
        {
            if (hasEdge(u, v))
                neighbours.push_back(v);
        }

        return neighbours;
    }
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
 *  Brute-force algorithm to find the Maximal Independent Set in a graph.
 *  This algorithm explores all possible subsets of vertices and checks if they form an independent set.
 */

int bruteForceRec(const Graph &g, int idx, std::vector<int> &current, std::vector<int> &best)
{
    int n = g.numVertices();

    if (idx == n)
    {
        if (isValidIndependentSet(g, current))
        {
            if (current.size() > best.size())
            {
                best = current; // Found a better independent set
            }
        }
        return best.size();
    }

    current.push_back(idx);
    bruteForceRec(g, idx + 1, current, best);
    current.pop_back();

    bruteForceRec(g, idx + 1, current, best);

    return best.size();
}

int bruteForce(const Graph &g, std::vector<int> &independentSet) {
    std::vector<int> current;
    bruteForceRec(g, 0, current, independentSet);
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
    int size = bruteForce(*graph, independentSet);
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