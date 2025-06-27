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

public:
  Graph(int numVertices) : n(numVertices)
  {
    adjMatrix = new bool *[n];
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
    }
  }

  void deleteEdge(int u, int v)
  {
    if (u != v && u >= 0 && v >= 0 && u < n && v < n)
    {
      adjMatrix[u][v] = adjMatrix[v][u] = false;
    }
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
  const bool* row(int u) const { return adjMatrix[u]; }
};

// ALGORITHM (START)
/**
 *  Algo v1 + pruning of partial solutions that can't become valid independent sets.
 */

void backtrackRec(const Graph &g, int idx, std::vector<int> &current, std::vector<int> &best)
{
    int n = g.numVertices();

    // if the last node has been considered return the better solution between the current and the best yet
    if (idx == n)
    {
        if (current.size() > best.size())
        {
            best = current;
        }
        return;
    }

    // check for conflicts of node idx with the current set
    bool canInclude = true;
    for (int v : current)
    {
        if (g.hasEdge(idx, v)) // conflict
        {
            canInclude = false;
            break;
        }
    }

    // include node idx
    if (canInclude)
    {
        current.push_back(idx);
        backtrackRec(g, idx + 1, current, best);
        current.pop_back();
    }

    // exclude node idx
    backtrackRec(g, idx + 1, current, best);
}

int backtrack(const Graph &g, std::vector<int> &independentSet)
{
    std::vector<int> current;
    backtrackRec(g, 0, current, independentSet);
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
    int size = backtrack(*graph, independentSet);
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