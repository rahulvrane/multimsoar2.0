# MultiMSOAR 2.0 - Code Structure & Threading Analysis

## FILE COMPOSITION BREAKDOWN

### 1. MultiMSOARSoftware.cpp (332 lines) - Main Entry Point
**Responsibility**: I/O handling, algorithm orchestration, result aggregation

Key sections:
```
Lines 15-37:    Global variable declarations (THREADING HAZARD)
Lines 39-49:    DFS() function - single linkage clustering
Lines 52-155:   Partition() function - Hungarian algorithm application
Lines 157-197:  TreeLabeling() function - algorithm selection & execution
Lines 199-217:  printGeneInfo() function - output generation
Lines 219-332:  main() function - orchestration & I/O
```

**Global State Issues** (Lines 15-37):
```cpp
map<string,int> visited;           // ← Used in DFS, NOT thread-safe
vector<string> group;              // ← Cleared per family, thread-local needed
set<string> AllGeneBirth;          // ← Needs mutex protection
set<string> AllGeneDuplication;    // ← Needs mutex protection
map<int, int> AllGeneLoss;         // ← Needs mutex protection
```

**Main Processing Loop** (Lines 309-326):
```cpp
// FAMILY PROCESSING - HIGHLY PARALLELIZABLE
for(map<int,set<string> >::iterator it=RealFamily.begin(); 
    it!=RealFamily.end(); it++)  // Iterate ~100+ families
{
    AllTrees.clear();              // Per-family clearing
    AllTreeGeneName.clear();       // Per-family clearing
    
    for(set<string>::iterator j=(it->second).begin(); 
        j!=(it->second).end(); j++)
    {
        if(visited.count(*j)==0)   // ← RACE CONDITION: visited[] must be thread-local
        {
            group.clear();         // ← RACE CONDITION: group[] must be thread-local
            DFS(*j);               // ← Process component
            Partition();           // ← Apply Hungarian algorithm
        }
    }
    TreeLabeling(orthoGroupOut);   // ← Choose algorithm & output
}
```

---

### 2. Hungarian.h (381 lines) - Bipartite Matching
**Algorithm**: Maximum weighted bipartite matching using augmenting paths
**Complexity**: O(n³) where n = matrix dimension
**Current Use**: Called once per internal node merge in Partition()

Key operations:
```
Constructor (lines 40-153):
  1. Initialize labeling (lines 56-70)          O(n)
  2. Generate equality subgraph (lines 72-80)   O(n²)
  3. Find initial maximal matching (lines 82-100) O(n²)
  4. Main loop (lines 104-148):                 O(n³)
     - Find free vertex
     - Update labels
     - Find augmenting path
```

**Parallelization Potential**: 
- Multiple Hungarian instances per family can run in parallel (limited by Partition tree structure)
- Not worth paralyzing single Hungarian instance (algorithm already optimized)

---

### 3. TreeAnalysis.h (421 lines) - Gene Event Classification
**Responsibility**: Tree traversal, ortholog group identification, gene event detection

Key classes:
```
NODE (lines 11-19):
  - Binary tree node structure
  - Stores: value, group index, left/right pointers
  
GeneInfo (lines 21-27):
  - Gene annotation: name + tree node index
  
TreeAnalysis (lines 29-421):
  - Constructor: Receives labeled tree set
  - Main methods:
    * findOrthoGroup() - traverse tree to identify groups
    * printOrthoGroups() - output ortholog groups
    * printGeneInfo() - identify gene events (birth/duplication/loss)
```

**Critical Loop - printGeneInfo** (Lines 314-420):
```cpp
// LAYER ITERATION - MODERATELY PARALLELIZABLE
for(int i=0; i<N; i++)                    // For each tree layer
{
    vector<GeneInfo* > stack;
    int index=0;
    
    for(int j=0; j<S; j++)                // For each position in tree
    {
        // ... complex tree traversal logic ...
        // Identifies:
        // - Gene births (lines 357-364)
        // - Gene duplications (lines 356-368)
        // - Gene losses (lines 383, 413)
    }
}
```

**Parallelization Strategy**:
- Each layer (i) can be processed independently
- Results accumulated with mutex protection
- Expected speedup: 2-8x (depending on N value)

---

### 4. NodeCentric.h (376 lines) - Small Tree Labeling (N < 5)
**Algorithm**: Dynamic programming on tree nodes
**Complexity**: Exponential in tree depth (acceptable for N < 5)

Key structure:
```
Node class (lines 13-22):
  - map<long long,int> changes      // Cost for each state
  - map<long long,long long> leftV  // Optimal left subtree value
  - map<long long,long long> rightV // Optimal right subtree value
  
Constructor (lines 223-374):
  1. Build tree from input (lines 234-242)
  2. For each internal node (lines 244-342):
     a. Compute all combinations (lines 268, 276)
     b. Merge left/right results (lines 281-335)
  3. Find minimum cost solution (lines 345-356)
  4. Trace back optimal labeling (lines 361-372)
```

**Performance Note**:
- Exponential in tree depth due to state space (lines 99-125: vp[i] combinations)
- Slower than TreeCentric but exact optimal solution
- Limited to N < 5 to keep runtime reasonable

---

### 5. TreeCentric.h (280 lines) - Large Tree Labeling (N >= 5)
**Algorithm**: Dynamic programming on entire trees
**Complexity**: O(S × T) where S = tree count, T = internal nodes

Key structure:
```
Tree class (lines 12-19):
  - map<long long,int> AccCost      // Accumulated cost
  - map<long long,int> CurLabel     // Current tree labeling
  - map<long long,long long> PreValue // Previous state
  
Constructor (lines 220-277):
  1. Initialize with first tree (line 231)
  2. For each tree (lines 235-239):
     a. Validate internal labeling (lines 190, 80-83)
     b. Update accumulated cost (lines 192-216)
     c. Check constraints:
        - Zero_One_Constraint (lines 145-185)
        - One_Oh_One_Constraint (lines 98-141)
  3. Trace back optimal solution (lines 259-266)
```

**Parallelization Potential**:
- UpdateCurrentTree() calls for each tree could run in parallel
- Each tree processed independently, then synchronized
- Expected speedup: 1.5-2x

---

## DATA FLOW DIAGRAM

```
┌─────────────────────────────────────────────────────────────┐
│                    INPUT LOADING (SEQUENTIAL)               │
├─────────────────────────────────────────────────────────────┤
│ 1. speciesTree file → Parse → Internal representation       │
│ 2. Si_Sj files     → Load → Build adjacency/edges maps      │
│ 3. geneFamily file → Parse → RealFamily map                 │
└────────────────────────┬────────────────────────────────────┘
                         │
                         ↓
┌─────────────────────────────────────────────────────────────┐
│              FAMILY PROCESSING (PARALLELIZABLE)             │
│  *** HIGHEST PRIORITY FOR MULTI-THREADING ***              │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│  [THREAD 1]          [THREAD 2]          [THREAD N]        │
│  ┌──────────────┐   ┌──────────────┐   ┌──────────────┐  │
│  │ Family 0     │   │ Family 1     │   │ Family 99    │  │
│  ├──────────────┤   ├──────────────┤   ├──────────────┤  │
│  │ DFS Cluster  │   │ DFS Cluster  │   │ DFS Cluster  │  │
│  ├──────────────┤   ├──────────────┤   ├──────────────┤  │
│  │ Partition    │   │ Partition    │   │ Partition    │  │
│  │ (Hungarian)  │   │ (Hungarian)  │   │ (Hungarian)  │  │
│  ├──────────────┤   ├──────────────┤   ├──────────────┤  │
│  │ TreeLabeling │   │ TreeLabeling │   │ TreeLabeling │  │
│  │ Analysis     │   │ Analysis     │   │ Analysis     │  │
│  ├──────────────┤   ├──────────────┤   ├──────────────┤  │
│  │ Results_0    │   │ Results_1    │   │ Results_99   │  │
│  └──────────────┘   └──────────────┘   └──────────────┘  │
│         │                  │                  │            │
└─────────┼──────────────────┼──────────────────┼────────────┘
          │                  │                  │
          └──────────────────┼──────────────────┘
                             │
                    SYNCHRONIZATION POINT
                   (Mutex-protected aggregation)
                             │
                             ↓
┌─────────────────────────────────────────────────────────────┐
│             RESULT AGGREGATION (SEQUENTIAL)                │
├─────────────────────────────────────────────────────────────┤
│ AllGeneBirth += Results[i].GeneBirth  (for all i)           │
│ AllGeneDuplication += Results[i].GeneDuplication            │
│ AllGeneLoss += Results[i].GeneLoss                          │
└────────────────────────┬────────────────────────────────────┘
                         │
                         ↓
┌─────────────────────────────────────────────────────────────┐
│                    OUTPUT GENERATION                        │
├─────────────────────────────────────────────────────────────┤
│ GeneInfo file: gene birth, duplication, loss summary        │
│ OrthoGroup file: ortholog group membership                 │
└─────────────────────────────────────────────────────────────┘
```

---

## LOOP ANALYSIS (93 total loops)

### Parallelizable Loops

#### Loop Type 1: Family Iteration (PRIMARY TARGET)
```cpp
// Location: MultiMSOARSoftware.cpp lines 309-326
for(map<int,set<string>>::iterator it=RealFamily.begin(); 
    it!=RealFamily.end(); it++)
```
- **Iterations**: ~100 families
- **Cost per iteration**: 10-100ms (highly variable)
- **Load balance**: Unbalanced (some families much larger)
- **Data sharing**: NONE (families independent)
- **Parallelization**: ✓✓✓ IDEAL

#### Loop Type 2: Gene Position Iteration in Component
```cpp
// Location: MultiMSOARSoftware.cpp lines 317-323
for(set<string>::iterator j=(it->second).begin(); 
    j!=(it->second).end(); j++)
```
- **Iterations**: Variable (1-50+ per family)
- **Cost per iteration**: 1-10ms
- **Load balance**: Balanced within family
- **Data sharing**: NONE within family
- **Parallelization**: Possible with careful visited[] management

#### Loop Type 3: Layer Iteration (SECONDARY TARGET)
```cpp
// Location: TreeAnalysis.h lines 327-419
for(int i=0; i<N; i++)  // N = number of layers
```
- **Iterations**: Variable (1-20, depends on duplication extent)
- **Cost per iteration**: 1-50ms (complexity depends on tree structure)
- **Load balance**: Likely unbalanced
- **Data sharing**: YES - GeneBirth, GeneDuplication, GeneLoss need mutex
- **Parallelization**: ✓✓ GOOD (2-8x speedup)

#### Loop Type 4: Matrix Operations (INSIDE HUNGARIAN)
```cpp
// Location: Hungarian.h lines 73-80 (Equality subgraph)
for (int i=0; i<n; i++) {
  for (int j=0; j<n; j++) {
    if (Lx[i]+Ly[j] == weight[i][j])
      El[i][j] = 1;
  }
}
```
- **Iterations**: O(n²) where n ≤ ~50 (matrix dimension)
- **Cost per iteration**: Trivial (single comparison)
- **Load balance**: Perfect
- **Data sharing**: NONE (disjoint El[i][j])
- **Parallelization**: ✓ POSSIBLE (minor benefit, limited by call frequency)

### Non-Parallelizable Loops

#### Loop Type 5: DFS Traversal (DEPENDENT STATE)
```cpp
// Location: MultiMSOARSoftware.cpp lines 43-48
for(int i=0; i<adjacency[cur].size(); i++) {
  string next=adjacency[cur][i];
  if(visited.count(next)==0)
    DFS(next);  // Recursive
}
```
- **Issue**: visited[] map modified during traversal
- **Problem**: Dependency chain prevents parallelization
- **Alternative**: Could parallelize different connected components, but overhead high

#### Loop Type 6: Tree Construction (STACK-BASED)
```cpp
// Location: TreeAnalysis.h lines 83-103
for(int i=0; i<treeLabel.size(); i++) {
  // ... stack-based tree construction ...
  // Sequential order REQUIRED for correctness
}
```
- **Issue**: Stack operations require sequential order
- **Problem**: Each iteration depends on previous
- **Parallelization**: ✗ NOT POSSIBLE

---

## PERFORMANCE CHARACTERISTICS

### Time Complexity Analysis

| Component | Complexity | Variable | Current Impact |
|-----------|-----------|----------|-----------------|
| DFS clustering | O(V+E) | V=genes, E=pairs | ~5-10% |
| Partition (Hungarian) | O(N³) × K | N=genes/species K=merges | ~40-60% |
| TreeLabeling (NodeCentric) | Exp(depth) | depth~log(N) | ~10-20% (for N<5) |
| TreeLabeling (TreeCentric) | O(S×T) | S=trees T=nodes | ~10-20% (for N≥5) |
| TreeAnalysis | O(N × tree_nodes) | N=layers | ~10-15% |
| **Total Processing** | **O(F × (HUN + TREE + ANA))** | F=families | **100%** |

### Expected Single-Family Processing Time
```
Small family (~10 genes):    5-10ms
Medium family (~30 genes):   20-100ms
Large family (~100 genes):   200-500ms+

Total for 100 families: 1-5 minutes (sequential)
With 8-core parallelization: 7-40 seconds (8x speedup)
```

---

## MEMORY USAGE ANALYSIS

### Per-Thread Memory (Family Processing Context)

```
Static per thread:
  ├─ species map (read-only from main)      ~100 KB
  ├─ adjacency map (read-only from main)    ~200 KB
  ├─ edges map (read-only from main)        ~100 KB
  
Dynamic per family:
  ├─ AllTrees vector (N layers × 100 bytes) ~1-10 KB
  ├─ AllTreeGeneName (N × S genes)          ~10-50 KB
  ├─ visited map (per family)               ~10-100 KB
  ├─ group vector (per family)              ~1-10 KB
  ├─ Hungarian matrices (N×N integers)      ~20-200 KB
  ├─ Node pointers (TreeAnalysis trees)     ~10-50 KB
  └─ TOTAL per thread                       ~500 KB - 2 MB

16 threads: ~8-32 MB additional overhead
```

### Shared Memory
```
Input data (loaded once):
  ├─ speciesTree string                     <1 KB
  ├─ RealFamily map (100+ families)         ~50 KB
  ├─ species map (1000+ genes)              ~100 KB
  ├─ adjacency map (10K+ edges)             ~200 KB
  ├─ edges map (10K+ pairs)                 ~100 KB
  └─ TOTAL                                  ~450 KB
```

**Conclusion**: Memory overhead from threading is minimal (~20 MB for 16 threads)

---

## SYNCHRONIZATION REQUIREMENTS

### Critical Sections (Minimal)

```cpp
// Synchronization Point 1: Result Aggregation
// Location: After each family processing
{
  lock_guard<mutex> lock(result_mutex);
  AllGeneBirth.insert(thread_GeneBirth.begin(), thread_GeneBirth.end());
  AllGeneDuplication.insert(thread_GeneDuplication.begin(), thread_GeneDuplication.end());
  for(auto& [species, count] : thread_GeneLoss)
    AllGeneLoss[species] += count;
}
// Expected lock time: ~0.1-1ms per family (negligible)

// Synchronization Point 2: Output File Writing
// Location: TreeLabeling() output
{
  lock_guard<mutex> lock(file_mutex);
  ta.printOrthoGroups(orthoGroupOut);
}
// Alternative: Per-thread output buffering (better)
```

### Lock Contention Analysis

```
Assumption: 8 threads, 100 families
  ├─ Result aggregation calls: 100 (assuming 1 per family)
  ├─ Average lock wait time: <1ms per thread
  ├─ Total synchronization time: ~10ms
  └─ Overhead as % of total: 0.1-1% (negligible)
```

**Conclusion**: Synchronization overhead is negligible; thread pool design sound

---

## RESOURCE UTILIZATION ESTIMATES

### Current (Single-Threaded)
```
CPU Usage: 1 core @ 100%
           7 cores @ 0%
           
Memory: ~50 MB for input data
        ~500 MB for temporary structures (worst case)
        
I/O: Sequential file reading (limited to disk I/O bandwidth)
```

### With 8-Core Parallelization
```
CPU Usage: 8 cores @ ~90-95% (good load distribution)
           Slight imbalance due to variable family sizes
           
Memory: ~50 MB shared (input data)
        ~16 MB thread overhead (1 MB × 8 threads)
        ~4 GB temporary structures (worst case, 500MB × 8)
        
I/O: Could read Si_Sj files in parallel (2-4x improvement)
```

**Recommendation**: 8-core parallelization is sweet spot
- Diminishing returns beyond 8 cores
- Memory overhead becomes significant (>16 GB for 32 cores)
- Load imbalance increases with more cores

