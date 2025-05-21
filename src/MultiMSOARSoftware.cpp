#include <iostream>
#include <map>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <set>
#include <utility> // For std::pair
#include <tuple>   // For std::tuple
#include <mutex>   // For std::mutex
#include <thread>  // For std::thread
#include <algorithm> // For std::min, std::max
#include "Hungarian.h"
#include "TreeCentric.h"
#include "NodeCentric.h"
#include "TreeAnalysis.h"

using namespace std;

// --- Global Variables ---
// S: Number of species in the species tree. Loaded from command line argument.
int S;
// speciesTree: String representation of the species tree in postfix notation. Parsed from input file.
// This structure guides the layer merging process in Partition_modified.
string speciesTree; 
// maximumN: Tracks the maximum number of layers (gene copies across species) found in any processed gene family.
// Updated sequentially after all threads complete, based on results from TreeLabeling_modified.
int maximumN; 

// species: Maps gene names to an integer ID (0 to S-1) representing their species.
// Loaded from pairwise orthology files (e.g., "S0_S1"). Assumed to be read-only during parallel processing phase.
map<string,int> species; 
// adjacency: Represents the gene graph where nodes are genes and edges indicate potential orthology or similarity.
// Used by DFS to find connected components (which are then treated as gene families). Loaded from pairwise orthology files.
// Assumed to be read-only after initial loading.
map<string, vector<string> > adjacency; 
// edges: Stores weights (e.g., similarity scores) between pairs of genes, used by the Hungarian algorithm
// within Partition_modified to optimize layer matching. Loaded from pairwise orthology files.
// Assumed to be read-only during parallel processing.
map<pair<string,string>, double> edges; 


// visited: Tracks visited genes during the initial DFS-based component collection phase.
// Ensures each gene is processed as part of only one component, preventing redundant work. Modified sequentially.
map<string,int> visited; 

// AllTrees: Legacy global variable. Originally used to store tree strings for a single gene group
// when processing was entirely sequential. Its use is superseded by local variables in `Partition_modified`
// and thread-specific storage (`TaskResult` in `main`) in the refactored, multithreaded design.
vector<string> AllTrees; 
// AllTreeGeneName: Legacy global variable, similar to AllTrees, for storing gene names corresponding to nodes.
// Superseded by local variables and thread-specific storage in the refactored design.
vector<vector<string> > AllTreeGeneName; 

// Global accumulators for results across all processed gene families.
// These are populated sequentially in Phase 3 of `main`, after all threads have completed their work and
// their individual results (from `TaskResult` structs) are aggregated.
// AllGeneBirth: Set of unique gene names identified as birth events across all processed families.
set<string> AllGeneBirth; 
// AllGeneDuplication: Set of unique gene names identified as duplication events across all processed families.
set<string> AllGeneDuplication; 
// AllGeneLoss: Map storing cumulative counts of gene loss events per species ID (key) across all processed families.
map<int, int> AllGeneLoss; 

// group: Global vector temporarily storing genes of the current connected component found by DFS.
// Modified sequentially during the component collection phase (Phase 1 in `main`).
vector<string> group; 

// --- Mutexes for Shared Resources ---
// ortho_out_mutex: Ensures thread-safe write access to the shared ortholog group output file stream (`orthoGroupOut`).
// Each thread created in `main` for processing a gene group will lock this mutex before calling
// `ta.printOrthoGroups()` (via `TreeLabeling_modified`) to prevent interleaved writes.
std::mutex ortho_out_mutex; 
// results_aggregation_mutex: Declared for potential future use if the final result aggregation phase (Phase 3 in `main`)
// were to be parallelized. Currently, aggregation is sequential, so this mutex is not actively used for that purpose.
// It serves as a placeholder for if direct updates to global accumulators by threads were implemented.
std::mutex results_aggregation_mutex; 

// Performs Depth First Search (DFS) to find connected components of genes based on orthology data (adjacency map).
// A connected component is treated as a gene family (or part of one) to be processed by the subsequent steps.
// This function modifies the global `group` vector (populating it with genes of the current component)
// and the global `visited` map (marking genes as visited to prevent redundant processing of the same component).
// It is called sequentially in `main` during Phase 1 (Component Collection).
void DFS(string cur)
{
	group.push_back(cur);
	visited[cur]=1;
	for(int i=0; i<adjacency[cur].size(); i++)
	{
		string next=adjacency[cur][i];
		if(visited.count(next)==0)
			DFS(next);
	}
}

// Partition each group into N layers
void Partition() // Original function, kept for now
{
	vector<vector<vector<string> > > v(S);

	for(int i=0; i<group.size(); i++)
	{
		int sp=species[group[i]];
		vector<string> tmp;
		tmp.push_back(group[i]);
		v[sp].push_back(tmp);
	}
	
	// N: the number of layers
	int N_local=0; // Renamed N to N_local to avoid conflict if S is also a parameter name
	for(int i=0; i<S; i++) N_local=max(N_local, (int)(v[i].size()));

	//cout<<"N="<<N_local<<endl;

	// Pad each species in a group with dummy vertices
	for(int i=0; i<S; i++)
	{
		vector<string> dummy;
		dummy.push_back("");

		for(int j=v[i].size(); j<N_local; j++) 
			v[i].push_back(dummy);
	}

	//////////////////////////////////
	
	vector<vector<vector<string> > > stack;

	int index=0;

	for(size_t i=0; i<speciesTree.size(); i++)
	{
		if(speciesTree[i]!='N')
		{
			stack.push_back(v[index++]);
		}
		else
		{
			vector<vector<string> > v2=stack.back(); stack.pop_back();
			vector<vector<string> > v1=stack.back(); stack.pop_back();

			vector<vector<int> > matrix(N_local, vector<int> (N_local) );

			// Calculate the added weight for an edge in Bipartite Graph
			for(int j=0; j<N_local; j++) for(int k=0; k<N_local; k++) matrix[j][k]=0;

			for(int j=0; j<N_local; j++) 
			{
				for(int k=0; k<N_local; k++)
				{
					for(size_t jj=0; jj<v1[j].size(); jj++)
					{
						const string& gene1 = v1[j][jj];
						if (gene1 == "") continue; 

						for(size_t kk=0; kk<v2[k].size(); kk++)
						{
							const string& gene2 = v2[k][kk];
							if (gene2 == "") continue; 

							auto edge_iter = edges.find(make_pair(gene1, gene2));
							if(edge_iter != edges.end())
							{
								matrix[j][k] += (int)edge_iter->second;
							}
						}
					}
				}
			}
			Hungarian H(matrix);
			for(int j=0; j<N_local; j++)
			{
				int p=H.matchingX[j];
				for(size_t k=0; k<v2[p].size(); k++)
					v1[j].push_back(v2[p][k]);
			}
			stack.push_back(v1);
		}
	}
	
	vector<string> local_AllTrees(N_local, string(S,'0'));
	vector<vector<string> > local_AllTreeGeneName(N_local, vector<string> (S));

	for(int i=0; i<N_local; i++) for(int j=0; j<S; j++) local_AllTreeGeneName[i][j]="";

	for(int i=0; i<N_local; i++)
	{
		for(size_t j=0; j<stack[0][i].size(); j++)
		{
			if(stack[0][i][j]!="") 
			{
				local_AllTrees[i][species[stack[0][i][j]]]='1';
				local_AllTreeGeneName[i][species[stack[0][i][j]]]=stack[0][i][j];
			}
		}
		for(size_t j=0; j<speciesTree.size(); j++) if(speciesTree[j]=='N')
			local_AllTrees[i]=local_AllTrees[i].substr(0, j)+"N"+local_AllTrees[i].substr(j);
	}

	for(int i=0; i<N_local; i++) AllTrees.push_back(local_AllTrees[i]);
	for(int i=0; i<N_local; i++) AllTreeGeneName.push_back(local_AllTreeGeneName[i]);
}

// Refactored Partition function. Operates on local data passed as parameters.
// This function takes a group of connected genes (current_group) and, using species information
// (species_map, species_tree_str, S_val) and gene similarity scores (edges_map),
// partitions the genes into N_for_partition layers. This layering aims to group orthologous genes
// across different species. The Hungarian algorithm is used to optimize layer merging during
// traversal of the species tree.
// This function is designed to be thread-safe as it only reads from shared global data (passed as const ref)
// and operates on local variables or output parameters.
//
// Parameters:
// - current_group: A vector of gene names belonging to one connected component (a potential gene family).
// - species_map: A map from gene names to their respective species ID (integer). (Read-only access)
// - species_tree_str: The species tree represented in postfix string notation. (Read-only access)
// - S_val: The total number of species being considered. (Read-only access)
// - edges_map: A map storing precomputed orthology scores (or other similarity measures) between pairs of genes. (Read-only access)
// - N_for_partition (output parameter): This will be set to the number of layers determined for this gene group.
//
// Returns:
// - A std::pair containing two elements:
//   1. std::vector<std::string>: A vector of tree strings. Each string represents a layer,
//      encoded in a format suitable for TreeCentric/NodeCentric algorithms.
//   2. std::vector<std::vector<std::string>>: A 2D vector storing gene names. For each layer (tree string),
//      it maps species IDs to the gene name present in that species for that layer.
std::pair<std::vector<std::string>, std::vector<std::vector<std::string>>> 
Partition_modified(
    const std::vector<std::string>& current_group, 
    const std::map<std::string,int>& species_map, 
    const std::string& species_tree_str, 
    int S_val, 
    const std::map<std::pair<std::string,std::string>, double>& edges_map,
    int& N_for_partition // Output parameter
)
{
	vector<vector<vector<string> > > v(S_val); // v[species_id] stores layers of genes for that species

	// Distribute genes in current_group into layers based on their species
	for(size_t i=0; i<current_group.size(); i++)
	{
        auto it = species_map.find(current_group[i]);
        if (it == species_map.end()) {
            cerr << "Warning: Gene " << current_group[i] << " not found in species_map. Skipping in Partition_modified." << endl;
            continue;
        }
		int sp = it->second;
		vector<string> tmp;
		tmp.push_back(current_group[i]);
		v[sp].push_back(tmp); // Each gene initially forms its own "layer" within its species
	}
	
	// Determine N_for_partition: the maximum number of layers (genes) in any single species for this group.
	N_for_partition=0;
	for(int i=0; i<S_val; i++) N_for_partition=max(N_for_partition, (int)(v[i].size()));

	// Pad layers with dummy empty strings so all species have N_for_partition layers.
	for(int i=0; i<S_val; i++)
	{
		vector<string> dummy;
		dummy.push_back("");
		for(int j=v[i].size(); j<N_for_partition; j++) 
			v[i].push_back(dummy);
	}
	
	// Process the species tree to merge layers using Hungarian algorithm for optimal matching.
	vector<vector<vector<string> > > stack; // Stack for species tree traversal (postfix)
	int index=0; // Index for accessing species data in v

	for(size_t i=0; i<species_tree_str.size(); i++)
	{
		if(species_tree_str[i]!='N') // Leaf node in species tree (represents a species)
		{
			stack.push_back(v[index++]);
		}
		else // Internal node in species tree (represents a speciation event)
		{
			// Pop two children (species or subtrees) from stack
			vector<vector<string> > v2=stack.back(); stack.pop_back();
			vector<vector<string> > v1=stack.back(); stack.pop_back();

			// Create cost matrix for Hungarian algorithm
			vector<vector<int> > matrix(N_for_partition, vector<int> (N_for_partition, 0) );

			// Populate cost matrix based on orthology scores between genes in corresponding layers
			for(int j=0; j<N_for_partition; j++) 
			{
				for(int k=0; k<N_for_partition; k++)
				{
					for(size_t jj=0; jj<v1[j].size(); jj++)
					{
						const string& gene1 = v1[j][jj];
						if (gene1 == "") continue; 

						for(size_t kk=0; kk<v2[k].size(); kk++)
						{
							const string& gene2 = v2[k][kk];
							if (gene2 == "") continue; 

							auto edge_iter = edges_map.find(make_pair(gene1, gene2));
							if(edge_iter != edges_map.end())
							{
								matrix[j][k] += (int)edge_iter->second;
							}
						}
					}
				}
			}
			// Find maximum weight bipartite matching to pair up layers
			Hungarian H(matrix);
			// Merge matched layers
			for(int j=0; j<N_for_partition; j++)
			{
				int p=H.matchingX[j]; // Layer p from v2 is matched with layer j from v1
				for(size_t k=0; k<v2[p].size(); k++)
					v1[j].push_back(v2[p][k]); // Add genes from v2's matched layer to v1's layer
			}
			stack.push_back(v1); // Push merged result back to stack
		}
	}
	
	// Construct final tree strings and gene name mappings
	std::vector<std::string> local_AllTrees(N_for_partition, string(S_val,'0'));
	std::vector<std::vector<std::string>> local_AllTreeGeneName(N_for_partition, vector<string> (S_val));

	// Initialize gene names to empty strings
	for(int i=0; i<N_for_partition; i++) for(int j=0; j<S_val; j++) local_AllTreeGeneName[i][j]="";

	// Populate tree strings and gene names based on the final merged stack
	for(int i=0; i<N_for_partition; i++) // For each layer
	{
		for(size_t j=0; j<stack[0][i].size(); j++) // For each gene in the layer
		{
			if(stack[0][i][j]!="") // If not a dummy gene
			{
                auto it = species_map.find(stack[0][i][j]);
                if (it != species_map.end()) {
				    local_AllTrees[i][it->second]='1'; // Mark presence in tree string
				    local_AllTreeGeneName[i][it->second]=stack[0][i][j]; // Store gene name
                } else {
                    cerr << "Warning: Gene " << stack[0][i][j] << " from stack not found in species_map during final tree construction." << endl;
                }
			}
		}
		// Insert 'N' for internal nodes into tree strings (consistent with TreeCentric/NodeCentric format)
		for(size_t j=0; j<species_tree_str.size(); j++) if(species_tree_str[j]=='N')
			local_AllTrees[i]=local_AllTrees[i].substr(0, j)+"N"+local_AllTrees[i].substr(j);
	}
    return std::make_pair(local_AllTrees, local_AllTreeGeneName); // Return results for this group
}


void TreeLabeling(ofstream& orthoGroupOut) // Original function, kept for now
{
	// Decide which labeling algorithm to use (NodeCentric or TreeCentric)
	vector<string> LabelResults;
	int totalSubstitutions;

	int N_local = AllTrees.size(); // Use local N
	maximumN = max(maximumN, N_local); // Still updates global maximumN

	if(N_local==0) return;

	if(N_local<5)
	{
		NodeCentric nc(AllTrees); // Uses global AllTrees
		totalSubstitutions=nc.totalSubstitutions;
		LabelResults=nc.optimalLabeling;
	}
	else
	{
		TreeCentric tc(AllTrees); // Uses global AllTrees
		totalSubstitutions=tc.totalSubstitutions;
		LabelResults=tc.optimalLabeling;
	}

	TreeAnalysis ta(speciesTree, LabelResults, AllTreeGeneName); // Uses global speciesTree, AllTreeGeneName
	ta.printOrthoGroups(orthoGroupOut);
	ta.printGeneInfo();
	for(size_t i=0; i<ta.GeneBirth.size(); i++) AllGeneBirth.insert(ta.GeneBirth[i]); // Updates global
	for(size_t i=0; i<ta.GeneDuplication.size(); i++) AllGeneDuplication.insert(ta.GeneDuplication[i]); // Updates global
	for(size_t i=0; i<ta.GeneLoss.size(); i++) AllGeneLoss[ta.GeneLoss[i]]++; // Updates global
}

// Refactored TreeLabeling function. Operates on local data passed as parameters.
// This function takes the tree strings and gene names (output from Partition_modified for a specific gene group),
// determines the optimal labeling of internal nodes using either NodeCentric or TreeCentric algorithms,
// and then performs tree analysis to identify gene birth, duplication, and loss events.
// Ortholog groups are written to a shared output stream, protected by a mutex for thread safety.
// This function is designed to be callable from multiple threads, each processing a different gene group.
//
// Parameters:
// - current_AllTrees: A vector of tree strings, where each string represents a layer and its gene presence/absence across species.
// - current_AllTreeGeneName: A 2D vector mapping each node in current_AllTrees to its corresponding gene name.
// - species_tree_str: The species tree represented in postfix string notation. (Read-only access)
// - ortho_out_stream: An output file stream to which ortholog group information is written. This stream is shared
//                     across threads, so access is controlled by ortho_out_mutex.
// - ortho_out_mutex: A mutex to ensure thread-safe writes to ortho_out_stream.
//
// Returns:
// - A std::tuple containing:
//   1. std::set<std::string>: A set of gene names identified as birth events within this gene group.
//   2. std::set<std::string>: A set of gene names identified as duplication events within this gene group.
//   3. std::map<int, int>: A map where keys are species IDs and values are counts of gene loss events in that species for this group.
//   4. int: The number of layers (N_val) processed for this gene group, which is current_AllTrees.size().
std::tuple<std::set<std::string>, std::set<std::string>, std::map<int, int>, int>
TreeLabeling_modified(
    const std::vector<std::string>& current_AllTrees, 
    const std::vector<std::vector<std::string>>& current_AllTreeGeneName, 
    const std::string& species_tree_str, 
    std::ofstream& ortho_out_stream, 
    std::mutex& ortho_out_mutex
)
{
    vector<string> LabelResults; // Stores the optimal labeling strings from NodeCentric/TreeCentric.
    int totalSubstitutions;      // Stores the total substitutions for the optimal labeling.

    int N_val = current_AllTrees.size(); // Number of layers/trees for this group.

    if(N_val == 0) { // If no trees were generated (e.g., empty group), return empty results.
        return std::make_tuple(std::set<std::string>(), std::set<std::string>(), std::map<int, int>(), 0);
    }

    // Choose labeling algorithm based on N_val (number of layers).
    // NodeCentric is typically faster for smaller N.
    if(N_val < 5)
    {
        NodeCentric nc(current_AllTrees);
        totalSubstitutions = nc.totalSubstitutions; 
        LabelResults = nc.optimalLabeling;
    }
    else
    {
        TreeCentric tc(current_AllTrees);
        totalSubstitutions = tc.totalSubstitutions; 
        LabelResults = tc.optimalLabeling;
    }

    // Perform tree analysis using the optimal labeling results.
    TreeAnalysis ta(species_tree_str, LabelResults, current_AllTreeGeneName);
    
    // Write ortholog groups to the shared output file stream, protected by a mutex.
    {
        std::lock_guard<std::mutex> lock(ortho_out_mutex);
        ta.printOrthoGroups(ortho_out_stream);
    }
    
    // This typically prints to std::cout or a similar general output, not a shared file per se.
    // If it were writing to a shared resource specific to this call that needed aggregation,
    // it would also need protection or a different handling strategy.
    ta.printGeneInfo(); 

    // Collect gene birth, duplication, and loss events for this group.
    std::set<std::string> local_GeneBirth;
    std::set<std::string> local_GeneDuplication;
    std::map<int, int> local_GeneLoss;

    for(size_t i=0; i<ta.GeneBirth.size(); i++) local_GeneBirth.insert(ta.GeneBirth[i]);
    for(size_t i=0; i<ta.GeneDuplication.size(); i++) local_GeneDuplication.insert(ta.GeneDuplication[i]);
    for(size_t i=0; i<ta.GeneLoss.size(); i++) local_GeneLoss[ta.GeneLoss[i]]++;
    
    // Return collected events and N_val.
    return std::make_tuple(local_GeneBirth, local_GeneDuplication, local_GeneLoss, N_val);
}


void printGeneInfo(char* filename)
{
	ofstream outfile(filename);
	outfile<<"Gene birth: ";
	for(set<string>::iterator it=AllGeneBirth.begin(); it!=AllGeneBirth.end(); it++)
		outfile<<*it<<"\t";
	outfile<<endl;

	outfile<<"Gene duplication: ";
	for(set<string>::iterator it=AllGeneDuplication.begin(); it!=AllGeneDuplication.end(); it++)
		outfile<<*it<<"\t";
	outfile<<endl;

	outfile<<"Gene loss: ";
	for(map<int,int>::iterator it=AllGeneLoss.begin(); it!=AllGeneLoss.end(); it++)
		outfile<<"Species"<<it->first<<"\t"<<it->second<<"\t";
	outfile<<endl;
	outfile.close();
}

int main(int argc, char** argv)
{
	if(argc!=6)
	{
		cout<<"Usage: MultiMSOAR2.0 <#species> <speciesTree> <GeneFamily> <-o GeneInfo> <-o OrthoGroups>"<<endl;
		exit(1);
	}

	ifstream infile(argv[2]);
	string tmpSpeciesTree;
	getline(infile, tmpSpeciesTree);
	infile.close();
	speciesTree="";
	string tmpL="";
	for(int i=0; i<tmpSpeciesTree.length(); i++)
	{
		if(tmpSpeciesTree[i]==',')
		{
			if(tmpL!="") { speciesTree+='1'; tmpL=""; }
		}
		else if(tmpSpeciesTree[i]==')')
		{
			if(tmpL!="")
			{
				speciesTree+='1';
			}
			speciesTree+='N';
			tmpL="";
			i++;
			while(i<tmpSpeciesTree.length() and tmpSpeciesTree[i]!=')' and tmpSpeciesTree[i]!=',')
				i++;
			i--;
		}
		else if(tmpSpeciesTree[i]==' ') continue;
		else
			tmpL+=tmpSpeciesTree[i];
	}
	

	S=atoi(argv[1]);
	//S=(speciesTree.size()+1)/2;

	int level=0, level2=1;
	for(int i=0; i<S; i++) for(int j=i+1; j<S; j++)
	{
		stringstream ss;
		ss<<"S"<<i<<"_S"<<j;
		ifstream file(ss.str().c_str());
		if( file.fail() )
		{
			cout<<"Cannot open file "<<ss.str()<<endl;
			exit(1);
		}

		string line;
		while(getline(file, line))
		{
			stringstream ss;
			ss<<line;
			string gene1, gene2;
			double score;
			ss>>gene1>>gene2>>score;

			species[gene1]=i;
			species[gene2]=j;

			adjacency[gene1].push_back(gene2);
			adjacency[gene2].push_back(gene1);

			edges[make_pair(gene1,gene2)]=score;
			edges[make_pair(gene2,gene1)]=score;
		}
	}
	
	map<int, set<string> > RealFamily;
	ifstream infile2(argv[3]);
	string oneFamily;
	int indexFamily=0;
	while(getline(infile2, oneFamily))
	{
		stringstream ss(oneFamily);
		string name;
		while(ss>>name) RealFamily[indexFamily].insert(name);
		indexFamily++;
	}
	infile2.close();


	// Sequential Component Collection
	std::vector<std::vector<std::string>> all_distinct_groups;
	// visited and group are still global and used by DFS as before.
	// This part remains sequential.
	for (auto it_family = RealFamily.begin(); it_family != RealFamily.end(); ++it_family) { // Using iterators for broader C++ version compatibility
		const std::set<std::string>& gene_set = it_family->second;
		for (const std::string& gene_name : gene_set) {
			if (visited.count(gene_name) == 0) {
				group.clear(); 
				DFS(gene_name); 
				if (!group.empty()) {
					all_distinct_groups.push_back(group); // Store a copy
				}
			}
		}
	}

	// --- Main Program Logic ---
	// Phase 1: Sequential Component Collection.
	// This phase identifies all distinct groups of connected genes (potential gene families)
	// from the input `RealFamily` data. It uses DFS, which modifies shared global state 
	// (`visited` map and `group` vector). Therefore, this part remains sequential to avoid
	// complex synchronization for DFS itself.
	std::vector<std::vector<std::string>> all_distinct_groups; // Stores all unique gene groups found.
	for (auto it_family = RealFamily.begin(); it_family != RealFamily.end(); ++it_family) {
		const std::set<std::string>& gene_set = it_family->second;
		for (const std::string& gene_name : gene_set) {
			if (visited.count(gene_name) == 0) { // If gene hasn't been visited yet
				group.clear();                   // Clear global group for this new component
				DFS(gene_name);                  // Perform DFS to find all genes in this component
				if (!group.empty()) {
					all_distinct_groups.push_back(group); // Store a copy of the found group
				}
			}
		}
	}

	// Phase 2: Parallel Processing of Components.
	// Each distinct gene group collected in Phase 1 is now processed in a separate thread.
	// This allows for concurrent execution of the computationally intensive Partition and TreeLabeling steps.
    std::vector<std::thread> worker_threads; // Vector to store thread objects.
    if (!all_distinct_groups.empty()) { 
        worker_threads.reserve(all_distinct_groups.size()); // Reserve space for efficiency.
    }

    // TaskResult struct to hold results from each thread's execution.
    // This allows each thread to work independently and store its findings without direct contention
    // on global result accumulators during the parallel phase.
    struct TaskResult {
        std::set<std::string> birth_events;
        std::set<std::string> duplication_events;
        std::map<int, int> loss_events;
        int n_val_for_max; // Number of layers from this task, for updating global maximumN.
        bool success = true; // Flag to indicate if the task (thread processing) completed successfully.
    };
    // Vector to store results from each task/thread. Sized to match the number of groups.
    std::vector<TaskResult> per_task_results(all_distinct_groups.size());

    ofstream orthoGroupOut(argv[5]); // Ortholog group output file, shared among threads. Opened once.

    // Launch one thread per distinct gene group.
    for (size_t i = 0; i < all_distinct_groups.size(); ++i) {
        worker_threads.emplace_back([/* Lambda Captures: */
                                     // Data specific to this task:
                                     &all_distinct_groups, i, // The vector of all groups (read-only for this thread's group) and current group index.
                                     // Read-only global data (captured by reference for efficiency, const correctness assumed by usage):
                                     &species, &speciesTree, S_param = S, &edges, 
                                     // Shared resources (must be handled carefully, e.g., with mutexes):
                                     &orthoGroupOut, &ortho_out_mutex, // Output stream and its mutex.
                                     // Output location for this thread's results:
                                     &per_task_results                 // Reference to vector where this thread stores its results.
                                     ] {
            try { // Add try-catch block for exception safety within the thread.
                const auto& current_component_group = all_distinct_groups[i];

                if (current_component_group.empty()) { // Should not happen if DFS ensures non-empty groups
                    if (i < per_task_results.size()) per_task_results[i].success = false;
                    return; // Skip empty groups.
                }

                int N_val_from_partition = 0; 
                
                // Step 1: Partition the current gene group. This is computationally intensive.
                std::pair<std::vector<std::string>, std::vector<std::vector<std::string>>> partition_output =
                    Partition_modified(current_component_group, species, speciesTree, S_param, edges, N_val_from_partition);

                // Step 2: Perform tree labeling and analysis. Also computationally intensive.
                std::tuple<std::set<std::string>, std::set<std::string>, std::map<int, int>, int> labeling_output =
                    TreeLabeling_modified(partition_output.first, partition_output.second, speciesTree, orthoGroupOut, ortho_out_mutex);

                // Store results from this thread's task into its designated slot in per_task_results.
                if (i < per_task_results.size()) { // Bounds check for safety.
                    per_task_results[i].birth_events = std::get<0>(labeling_output);
                    per_task_results[i].duplication_events = std::get<1>(labeling_output);
                    per_task_results[i].loss_events = std::get<2>(labeling_output);
                    per_task_results[i].n_val_for_max = std::get<3>(labeling_output);
                    per_task_results[i].success = true;
                }

            } catch (const std::exception& e) { // Catch standard exceptions.
                std::cerr << "Thread for group " << i << " (first gene: " 
                          << (all_distinct_groups[i].empty() ? "N/A" : all_distinct_groups[i][0]) 
                          << ") caught an exception: " << e.what() << std::endl;
                if (i < per_task_results.size()) per_task_results[i].success = false; // Mark task as failed
            } catch (...) { // Catch any other type of exception.
                std::cerr << "Thread for group " << i << " (first gene: " 
                          << (all_distinct_groups[i].empty() ? "N/A" : all_distinct_groups[i][0]) 
                          << ") caught an unknown exception." << std::endl;
                if (i < per_task_results.size()) per_task_results[i].success = false; // Mark task as failed
            }
        });
    }

    // Wait for all worker threads to complete their execution.
    for (std::thread& t : worker_threads) {
        if (t.joinable()) {
            t.join();
        }
    }
    
    orthoGroupOut.close(); // Close the shared output file stream after all threads are done.

    // Phase 3: Sequential Result Aggregation.
    // After all threads have finished, collect results from per_task_results.
    // This aggregation is done sequentially to avoid needing the results_aggregation_mutex
    // for the global accumulators (AllGeneBirth, etc.), simplifying logic for this step.
    maximumN = 0; // Reset global maximumN before aggregation.
    AllGeneBirth.clear(); 
    AllGeneDuplication.clear(); 
    AllGeneLoss.clear();

    for (size_t k = 0; k < per_task_results.size(); ++k) {
        if (per_task_results[k].success) { // Only aggregate if the task was successful.
            AllGeneBirth.insert(per_task_results[k].birth_events.begin(), per_task_results[k].birth_events.end());
            AllGeneDuplication.insert(per_task_results[k].duplication_events.begin(), per_task_results[k].duplication_events.end());
            for(const auto& pair_val : per_task_results[k].loss_events) {
                AllGeneLoss[pair_val.first] += pair_val.second;
            }
            maximumN = std::max(maximumN, per_task_results[k].n_val_for_max);
        }
    }
	
	printGeneInfo(argv[4]);

	return 0;
}
