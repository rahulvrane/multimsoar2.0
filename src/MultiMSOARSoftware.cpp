#include <iostream>
#include <map>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <set>
#include <thread>
#include <mutex>
#include <future>
#include <chrono>
#include <algorithm>
#include "Hungarian.h"
#include "TreeCentric.h"
#include "NodeCentric.h"
#include "TreeAnalysis.h"
#include "ThreadPool.h"

using namespace std;

// Species tree and the number of species
int S;
string speciesTree;
int maximumN;

// Vertex and Adjacency list
map<string,int> species;
map<string, vector<string> > adjacency;
map<pair<string,string>, double> edges;


// Global results (thread-safe access via ResultAggregator)
set<string> AllGeneBirth;
set<string> AllGeneDuplication;
map<int, int> AllGeneLoss;

// Mutex for output file writing (used in parallel processing)
mutex output_mutex;

// Thread-local DFS - accepts local visited and group to avoid race conditions
void DFS_Local(string cur,
               map<string,int>& visited_local,
               vector<string>& group_local,
               const map<string, vector<string> >& adjacency)
{
	group_local.push_back(cur);
	visited_local[cur]=1;
	for(int i=0; i<adjacency.at(cur).size(); i++)
	{
		string next=adjacency.at(cur)[i];
		if(visited_local.count(next)==0)
			DFS_Local(next, visited_local, group_local, adjacency);
	}
}

// Partition each group into N layers (thread-safe version)
void Partition_Local(const vector<string>& group_local,
                     vector<string>& AllTrees_local,
                     vector<vector<string> >& AllTreeGeneName_local,
                     const map<string,int>& species,
                     const map<pair<string,string>, double>& edges,
                     const string& speciesTree,
                     int S)
{
	vector<vector<vector<string> > > v(S);

	for(int i=0; i<group_local.size(); i++)
	{
		int sp=species.at(group_local[i]);
		vector<string> tmp;
		tmp.push_back(group_local[i]);
		v[sp].push_back(tmp);
	}
	
	// N: the number of layers
	int N=0;
	for(int i=0; i<S; i++) N=max(N, (int)(v[i].size()));

	//cout<<"N="<<N<<endl;

	// Pad each species in a group with dummy vertices
	for(int i=0; i<S; i++)
	{
		vector<string> dummy;
		dummy.push_back("");

		for(int j=v[i].size(); j<N; j++) 
			v[i].push_back(dummy);
	}

	//////////////////////////////////
	
	vector<vector<vector<string> > > stack;

	int index=0;

	for(int i=0; i<speciesTree.size(); i++)
	{
		if(speciesTree[i]!='N')
		{
			stack.push_back(v[index++]);
		}
		else
		{
			vector<vector<string> > v2=stack.back(); stack.pop_back();
			vector<vector<string> > v1=stack.back(); stack.pop_back();

			vector<vector<int> > matrix(N, vector<int> (N) );

			// Calculate the added weight for an edge in Bipartite Graph
			for(int j=0; j<N; j++) for(int k=0; k<N; k++) matrix[j][k]=0;

			for(int j=0; j<N; j++) for(int k=0; k<N; k++)
			{
				for(int jj=0; jj<v1[j].size(); jj++) for(int kk=0; kk<v2[k].size(); kk++)
				{
					string gene1=v1[j][jj];
					string gene2=v2[k][kk];
					if(gene1=="" or gene2=="") matrix[j][k]+=0;
					auto edge_key = make_pair(gene1, gene2);
					if(edges.count(edge_key)) matrix[j][k]+=(int)edges.at(edge_key);
				}
			}

			// Run the Hungarian maximum matching algorithm for weighted bipartite graph
			//cout<<"Running Hungarian ..."<<endl;
			Hungarian H(matrix);
			//cout<<"Hungarian Done."<<endl;

			// Merge the vertices after matching
			for(int j=0; j<N; j++)
			{
				int p=H.matchingX[j];
				for(int k=0; k<v2[p].size(); k++)
					v1[j].push_back(v2[p][k]);
			}

			stack.push_back(v1);
		}
	}
	
	// Cout the partition for each group
	//cout<<N<<" layers: "<<endl;
	vector<string> trees(N, string(S,'0'));
	vector<vector<string> > treeGeneName(N, vector<string> (S));

	for(int i=0; i<N; i++) for(int j=0; j<S; j++) treeGeneName[i][j]="";

	for(int i=0; i<N; i++)
	{
		for(int j=0; j<stack[0][i].size(); j++)
		{
			if(stack[0][i][j]!="")
			{
				int sp = species.at(stack[0][i][j]);
				trees[i][sp]='1';
				treeGeneName[i][sp]=stack[0][i][j];
			}
		}
		for(int j=0; j<speciesTree.size(); j++) if(speciesTree[j]=='N')
			trees[i]=trees[i].substr(0, j)+"N"+trees[i].substr(j);
		//cout<<trees[i]<<endl;
	}

	// Store the results in local storage (thread-safe)
	for(int i=0; i<N; i++) AllTrees_local.push_back(trees[i]);
	for(int i=0; i<N; i++) AllTreeGeneName_local.push_back(treeGeneName[i]);
}

// Thread-safe TreeLabeling - accumulates results locally
void TreeLabeling_Local(const vector<string>& AllTrees_local,
                        const vector<vector<string> >& AllTreeGeneName_local,
                        set<string>& GeneBirth_local,
                        set<string>& GeneDuplication_local,
                        map<int, int>& GeneLoss_local,
                        stringstream& orthoGroupBuffer,
                        const string& speciesTree)
{
	// Decide which labeling algorithm to use (NodeCentric or TreeCentric)
	vector<string> LabelResults;
	int totalSubstitutions;

	int N=AllTrees_local.size();

	if(N==0) return;

	if(N<5)
	{
		NodeCentric nc(AllTrees_local);
		totalSubstitutions=nc.totalSubstitutions;
		LabelResults=nc.optimalLabeling;
	}
	else
	{
		TreeCentric tc(AllTrees_local);
		totalSubstitutions=tc.totalSubstitutions;
		LabelResults=tc.optimalLabeling;
	}

	//cout<<"Optimal Labeling: "<<totalSubstitutions<<endl;
	//for(int i=0; i<N; i++)
	//	cout<<LabelResults[i]<<endl;

	// Tree Analysis
	//cout<<"Tree Analysis: "<<endl;

	TreeAnalysis ta(speciesTree, LabelResults, AllTreeGeneName_local);
	//ta.printAnalysis();
	//ta.printDetailedAnalysis();

	// Write to buffer instead of file (for thread safety)
	ta.printOrthoGroups_Buffer(orthoGroupBuffer);

	ta.printGeneInfo();

	// Accumulate results locally
	for(int i=0; i<ta.GeneBirth.size(); i++) GeneBirth_local.insert(ta.GeneBirth[i]);
	for(int i=0; i<ta.GeneDuplication.size(); i++) GeneDuplication_local.insert(ta.GeneDuplication[i]);
	for(int i=0; i<ta.GeneLoss.size(); i++) GeneLoss_local[ta.GeneLoss[i]]++;
	//cout<<"********************************************************"<<endl;
}

/**
 * Process a single gene family (thread-safe worker function)
 * This function is called by the thread pool for parallel processing
 */
void processFamilyTask(int family_id,
                       const set<string>& family_genes,
                       const map<string,int>& species,
                       const map<string, vector<string> >& adjacency,
                       const map<pair<string,string>, double>& edges,
                       const string& speciesTree,
                       int S,
                       ResultAggregator& aggregator,
                       ofstream& orthoGroupOut)
{
	// Thread-local state
	map<string,int> visited_local;
	set<string> GeneBirth_local;
	set<string> GeneDuplication_local;
	map<int, int> GeneLoss_local;
	stringstream orthoGroupBuffer;

	// Process each gene in the family
	for(set<string>::const_iterator j=family_genes.begin(); j!=family_genes.end(); j++)
	{
		if(visited_local.count(*j)==0)
		{
			vector<string> group_local;
			vector<string> AllTrees_local;
			vector<vector<string> > AllTreeGeneName_local;

			// Perform DFS to find connected components
			DFS_Local(*j, visited_local, group_local, adjacency);

			// Partition the group into layers
			Partition_Local(group_local, AllTrees_local, AllTreeGeneName_local,
			               species, edges, speciesTree, S);

			// Perform tree labeling and analysis
			TreeLabeling_Local(AllTrees_local, AllTreeGeneName_local,
			                  GeneBirth_local, GeneDuplication_local, GeneLoss_local,
			                  orthoGroupBuffer, speciesTree);
		}
	}

	// Aggregate results (thread-safe)
	aggregator.aggregate(GeneBirth_local, GeneDuplication_local, GeneLoss_local);

	// Write output buffer to file (thread-safe)
	{
		lock_guard<mutex> lock(output_mutex);
		orthoGroupOut << orthoGroupBuffer.str();
	}
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

	cout << "MultiMSOAR 2.0 - Multi-threaded Edition" << endl;
	cout << "Hardware threads available: " << thread::hardware_concurrency() << endl;

	// ===========================================
	// TIER 2: PARALLEL FILE I/O
	// ===========================================
	cout << "Loading ortholog pair files in parallel..." << endl;
	auto io_start = chrono::high_resolution_clock::now();

	// Create thread pool for file I/O
	ThreadPool io_pool(min(8, (int)thread::hardware_concurrency()));

	// Vector to store file I/O tasks
	vector<future<IOTask>> io_futures;

	// Launch parallel file reads
	for(int i=0; i<S; i++) {
		for(int j=i+1; j<S; j++) {
			io_futures.push_back(io_pool.enqueue([i, j]() -> IOTask {
				IOTask task(i, j, "");
				stringstream ss;
				ss<<"S"<<i<<"_S"<<j;
				task.filename = ss.str();

				ifstream file(task.filename.c_str());
				if(file.fail()) {
					cerr<<"Cannot open file "<<task.filename<<endl;
					return task;
				}

				string line;
				while(getline(file, line)) {
					stringstream line_ss;
					line_ss<<line;
					string gene1, gene2;
					double score;
					line_ss>>gene1>>gene2>>score;

					task.ortholog_pairs.push_back(make_tuple(gene1, gene2, score));
				}
				return task;
			}));
		}
	}

	// Wait for all I/O tasks to complete and aggregate results
	for(auto& fut : io_futures) {
		IOTask task = fut.get();
		if(task.ortholog_pairs.empty()) {
			cout<<"Warning: No data loaded from file "<<task.filename<<endl;
			continue;
		}

		// Merge into global data structures
		for(const auto& pair : task.ortholog_pairs) {
			string gene1 = get<0>(pair);
			string gene2 = get<1>(pair);
			double score = get<2>(pair);

			species[gene1] = task.i;
			species[gene2] = task.j;

			adjacency[gene1].push_back(gene2);
			adjacency[gene2].push_back(gene1);

			edges[make_pair(gene1,gene2)] = score;
			edges[make_pair(gene2,gene1)] = score;
		}
	}

	auto io_end = chrono::high_resolution_clock::now();
	auto io_duration = chrono::duration_cast<chrono::milliseconds>(io_end - io_start);
	cout << "File I/O completed in " << io_duration.count() << " ms" << endl;
	
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


	// ===========================================
	// TIER 1: PARALLEL FAMILY PROCESSING
	// ===========================================
	cout << "Processing " << RealFamily.size() << " gene families in parallel..." << endl;
	auto family_start = chrono::high_resolution_clock::now();

	// Create result aggregator
	ResultAggregator aggregator(AllGeneBirth, AllGeneDuplication, AllGeneLoss);

	// Open output file
	ofstream orthoGroupOut(argv[5]);

	// Create thread pool for family processing
	ThreadPool family_pool(thread::hardware_concurrency());

	// Launch parallel family processing
	vector<future<void>> family_futures;

	// Create local references to avoid lambda capture warnings
	const map<string,int>& species_ref = species;
	const map<string, vector<string>>& adjacency_ref = adjacency;
	const map<pair<string,string>, double>& edges_ref = edges;
	const string& speciesTree_ref = speciesTree;
	const int S_val = S;

	for(map<int,set<string> >::iterator it=RealFamily.begin(); it!=RealFamily.end(); it++)
	{
		int family_id = it->first;
		set<string> family_genes = it->second;

		// Enqueue family processing task
		family_futures.push_back(family_pool.enqueue([family_id, family_genes,
		                                              &species_ref, &adjacency_ref,
		                                              &edges_ref, &speciesTree_ref,
		                                              S_val, &aggregator, &orthoGroupOut]() {
			processFamilyTask(family_id, family_genes, species_ref, adjacency_ref, edges_ref,
			                 speciesTree_ref, S_val, aggregator, orthoGroupOut);
		}));
	}

	// Wait for all families to complete
	for(auto& fut : family_futures) {
		fut.get();
	}

	orthoGroupOut.close();

	auto family_end = chrono::high_resolution_clock::now();
	auto family_duration = chrono::duration_cast<chrono::milliseconds>(family_end - family_start);
	cout << "Family processing completed in " << family_duration.count() << " ms" << endl;

	// Print summary statistics
	auto total_duration = chrono::duration_cast<chrono::milliseconds>(family_end - io_start);
	cout << "Total execution time: " << total_duration.count() << " ms" << endl;
	cout << "Gene birth events: " << AllGeneBirth.size() << endl;
	cout << "Gene duplication events: " << AllGeneDuplication.size() << endl;
	cout << "Gene loss events: " << AllGeneLoss.size() << endl;

	printGeneInfo(argv[4]);

	return 0;
}
