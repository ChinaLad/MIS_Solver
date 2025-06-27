 #include <bits/stdc++.h>
 #pragma GCC optimize("O3,unroll-loops")
 #pragma GCC target("sse4.2,popcnt")
 
 using namespace std;
 
 static constexpr int MAXN = 2048;     // raise if needed
 using B = std::bitset<MAXN>;          // handy alias
 
 class Graph {
 public:
     int           n;                  // #vertices
     long long     m = 0;              // #edges
     vector<B>     adj;                // (original)
     vector<B>     adjC;               // (complement)
     vector<int>   deg;                // degree(v)
 
     explicit Graph(int n_)
         : n(n_), adj(n_), adjC(n_), deg(n_, 0)
     {
         if (n > MAXN) throw runtime_error("MAXN too small for this graph");
         for (int i = 0; i < n; ++i) {
             adj[i].reset();
             adjC[i].reset().flip();          // all 1s
             adjC[i].reset(i);                // remove self-loop
         }
     }
 
     static Graph fromFile(const string& file)
     {
         ifstream in(file);
         if (!in) throw runtime_error("cannot open " + file);
         int N; in >> N;
         Graph g(N);
         for (int u, v; in >> u >> v;) g.addEdge(u, v);
         return g;
     }
 
     void addEdge(int u, int v)
     {
         if (u == v || u < 0 || v < 0 || u >= n || v >= n) return;
         if (adj[u].test(v)) return;          // already present
         adj [u].set(v);  adj [v].set(u);
         adjC[u].reset(v); adjC[v].reset(u);
         ++deg[u]; ++deg[v]; ++m;
     }
 
     /* handy row accessors */
     const B& nbr (int v) const { return adj [v]; }
     const B& cnbr(int v) const { return adjC[v]; }
 };
 

//Bron–Kerbosch with Tomita pivot  
 static int choosePivot(const Graph& G, const B& P, const B& X)
 {
     B PX = P | X;
     int bestU = -1, bestCnt = -1;
     for (int u = PX._Find_first(); u < G.n; u = PX._Find_next(u)) {
         int c = (P & G.cnbr(u)).count();     // |P ∩ complement(u)|
         if (c > bestCnt) { bestCnt = c; bestU = u; }
     }
     return bestU;
 }
 
 static void BK_MIS(const Graph& G,
                    B R, B P, B X,
                    vector<int>& best)
 {
     if (P.none() && X.none()) {              // R is maximal
         if (R.count() > best.size()) {
             best.clear();
             for (int v = R._Find_first(); v < G.n; v = R._Find_next(v))
                 best.push_back(v);
         }
         return;
     }
     int u = choosePivot(G, P, X);
     B ext = P & ~G.cnbr(u);                  // Tomita: P \ complement(u)
 
     for (int v = ext._Find_first(); v < G.n; v = ext._Find_next(v)) {
         BK_MIS(G,
                R | B{}.set(v),
                P & G.cnbr(v),
                X & G.cnbr(v),
                best);
         P.reset(v);
         X.set(v);
     }
 }
 

   //Greedy clique-cover bound  (colouring of the complement)
 static int cliqueCoverUB(const Graph& G, B cand)
 {
     int colours = 0;
     while (cand.any()) {
         B clique;
         int v = cand._Find_first();
         clique.set(v);
         cand.reset(v);
 
         for (int u = cand._Find_first(); u < G.n; u = cand._Find_next(u)) {
             /* add u if it is adjacent to EVERY vertex already in clique
              * ->   no non-edge in the original -> (complement(u) ∩ clique) = ∅ */
             if ((G.cnbr(u) & clique).none()) {
                 clique.set(u);
                 cand.reset(u);
             }
         }
         ++colours;                           // one more clique in the cover
     }
     return colours;
 }
 
 
 // Branch-and-Bound (sparse graphs)  –  exact
 static void BB_MIS(const Graph&  G,
                    B             cand,
                    B             chosen,
                    vector<int>&  best)
 {
     /* trivial & clique-cover bounds -------------------------------- */
     if (chosen.count() + cand.count()            <= best.size()) return;
     if (chosen.count() + cliqueCoverUB(G, cand)  <= best.size()) return;
 
     /* finished branch ---------------------------------------------- */
     if (cand.none()) {
         if (chosen.count() > best.size()) {
             best.clear();
             for (int v = chosen._Find_first(); v < G.n; v = chosen._Find_next(v))
                 best.push_back(v);
         }
         return;
     }
 
     // pick branching vertex = highest degree within cand
     int v = -1, bestDeg = -1;
     for (int u = cand._Find_first(); u < G.n; u = cand._Find_next(u)) {
         int d = (G.nbr(u) & cand).count();
         if (d > bestDeg) { bestDeg = d; v = u; }
     }
 
     //include v
     BB_MIS(G,
            cand & ~G.nbr(v) & ~B{}.set(v),
            chosen | B{}.set(v),
            best);
 
     // exclude v 
     cand.reset(v);
     BB_MIS(G, cand, chosen, best);
 }
 
//top-level driver
 int maxIndependentSet(const Graph& G, vector<int>& sol)
 {
     const double dens = (2.0 * G.m) / (G.n * (G.n - 1)); // edge density
     B cand; for (int i = 0; i < G.n; ++i) cand.set(i);
 
     vector<int> best;
     if (dens >= 0.6) {            // dense graph : clique in complement
         BK_MIS(G, {}, cand, {}, best);
     } else {                      //sparse graph : B&B
         BB_MIS(G, cand, {}, best);
     }
     sort(best.begin(), best.end());
     sol.swap(best);
     return static_cast<int>(sol.size());
 }
 
 int main(int argc, char* argv[])
 {
     if (argc < 2) {
         cerr << "Usage: " << argv[0] << "  <graph-file>\n";
         return 1;
     }
     Graph G = Graph::fromFile(argv[1]);
 
     vector<int> independentSet;
     auto t0 = chrono::high_resolution_clock::now();
     int  size = maxIndependentSet(G, independentSet);
     auto t1 = chrono::high_resolution_clock::now();
 
     cout << "MIS size: " << size << "\nVertices:";
     for (int v : independentSet) cout << ' ' << v;
     cout << "\nExecution time: "
          << chrono::duration<double>(t1 - t0).count()
          << " seconds\n";
     return 0;
 }
 