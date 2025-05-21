#include <iostream>
#include <algorithm>
#include <map>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <set>
#include <vector> // Added for std::vector
#include <algorithm> // Added for std::lexicographical_compare

using namespace std;

// Define Key Type and Comparator for TreeCentric algorithm.
// StateKey represents the binary labeling (0 or 1) of all S internal nodes of a tree.
// Using std::vector<bool> allows for S > 64, overcoming the limit of a long long bitmask.
// Each element in the vector corresponds to an internal node's label.
using StateKey = std::vector<bool>; 

// Custom comparator for StateKey to be used in std::map.
// Enables StateKey instances to be used as keys in ordered maps.
struct StateKeyCompare {
    bool operator()(const StateKey& a, const StateKey& b) const {
        // Lexicographical comparison is a standard way to compare sequences.
        return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
    }
};

class Tree
{
	public:
		map<StateKey, int, StateKeyCompare> AccCost;
		map<StateKey, StateKey, StateKeyCompare> CurLabel;
		map<StateKey, StateKey, StateKeyCompare> PreValue;
		Tree(){}
};


class TreeCentric
{
	public:

// Store the labeling results
int totalSubstitutions;
vector<string> optimalLabeling;

// N: the number of trees
// S: the number of internal nodes
int N, S;

// All the trees
vector<string> trees;

// Internal vector to store the labeling of the current tree
vector<int> w;

private:
    std::vector<Tree*> tree_layers; // Stores the Tree objects for each layer

public:
    // Constructor
    TreeCentric(vector<string> input); 
    // Destructor: Responsible for cleaning up dynamically allocated Tree objects
    // stored in the tree_layers vector to prevent memory leaks.
    ~TreeCentric();

// Given a tree, the valid labeling of the tree and its corresponding
// cost is stored in validLabeling
map<StateKey, int, StateKeyCompare> validLabeling;

// Recursively explores possible labelings for internal nodes of a single tree.
// 'tree': The current state of the tree string being processed (postfix notation).
// 'current_label': The StateKey being built, representing labels of internal nodes encountered so far.
// 'cost': The accumulated cost (number of substitutions) for the current_label.
void go2(string tree, StateKey current_label, int cost) // S is member variable TreeCentric::S
{
    if(tree.size()==1)
    {
        // Ensure label is of size S, padding if necessary (though ideally it should be built correctly)
        // This case implies that the label represents the full set of internal nodes.
        // The original code implies the label is built up one bit per 'N' encountered.
        // If current_label.size() != S, there's a mismatch in understanding or S needs to be passed/used.
        // For now, assume current_label is correctly built to size S by the time tree.size() == 1.
        if(validLabeling.count(current_label)==0 || validLabeling[current_label]>cost)
            validLabeling[current_label]=cost;
        return;
    }

    int N_pos=tree.find('N');
    if (N_pos == string::npos) { // Should not happen if tree.size() > 1 and logic is correct
        // Handle error or unexpected state
        return;
    }
    int left=tree[N_pos-2]-'0';
    int right=tree[N_pos-1]-'0';
    string prefix=tree.substr(0, N_pos-2);
    string suffix=tree.substr(N_pos+1);

    if(left==0 && right==0)
    {
        StateKey next_label = current_label;
        next_label.push_back(false);
        go2(prefix+string(1,'0')+suffix, next_label, cost);
    }
    else if(left==1 && right==1)
    {
        StateKey next_label = current_label;
        next_label.push_back(true);
        go2(prefix+string(1,'1')+suffix, next_label, cost);
    }
    else if((left==1 && right==0) || (left==0 && right==1))
    {
        StateKey next_label_1 = current_label;
        next_label_1.push_back(true);
        go2(prefix+string(1,'1')+suffix, next_label_1, cost+1);

        StateKey next_label_0 = current_label;
        next_label_0.push_back(false);
        go2(prefix+string(1,'2')+suffix, next_label_0, cost+1);
    }
    else // Handles cases like (2,0), (0,2), (2,1), (1,2), (2,2) etc.
         // This is the case where the resulting character in the tree string is '2'.
         // Original: go2(prefix+string(1,'2')+suffix, (label<<1)|0, cost+((left|right)&1));
    {
        StateKey next_label = current_label;
        next_label.push_back(false); // Appends 0 to the label, as per (label<<1)|0
        go2(prefix+string(1,'2')+suffix, next_label, cost+((left|right)&1));
    }
}

void Valid_Internal_Labeling(string tree)
{
    validLabeling.clear();
    StateKey initial_label; // Empty vector
    // The label should grow to size S. S is the number of 'N's in the original tree string.
    go2(tree, initial_label, 0);
    // cout<<"Valid Labeling: "<<validLabeling.size()<<endl;
    // After go2 completes, all keys in validLabeling should have size S.
    // We might need to verify this or ensure S is used correctly in go2 if padding is needed.
    // However, the recursive structure of go2 suggests it naturally creates labels of length S.
}

// Print the bit value of a StateKey
string printValue(const StateKey& p) // S is a member of TreeCentric
{
    string s = "";
    // S should be equal to p.size()
    // If p can have a different size than S, then the loop limit needs to be p.size()
    // Assuming p will always have size S when this is called.
    for(int i = 0; i < p.size(); ++i) // Iterate from the beginning of the vector
    {
        s += (p[i] ? '1' : '0');
    }
    return s;
}

// Check 1-0*-1 constraint between trees (InterTree 10*1 constraint)
// Parameter p is the StateKey representing the labeling of internal nodes.
bool One_Oh_One_Constraint(const StateKey& p_state) // p_state is the new name for the StateKey parameter
{
    // S is the number of internal nodes, p_state.size() should be S.
    // The original code reversed the bits from long long.
    // StateKey is already in the correct order [b0, b1, ..., bS-1]
    // where b0 is the label for the first 'N' encountered, etc.

    vector<int> stack;
    int index = 0; // Index for accessing elements of p_state
    for(int i=0; i<trees[0].size(); i++) // Iterate through the tree structure string
    {
        if(trees[0][i]=='0' || trees[0][i]=='1')
        {
            int tmp=0;
            for(int j=0; j<N; j++) // N is the number of trees
                tmp+=trees[j][i]-'0';

            if(tmp>0) stack.push_back(1); // Represents an aggregated '1'
            else stack.push_back(0);      // Represents an aggregated '0'
        }
        else // Encounter 'N', an internal node
        {
            if (index >= p_state.size()) {
                // This would indicate an error: more internal nodes in tree structure than in p_state
                // Or S is not correctly representing the number of internal nodes for p_state
                // For now, assume p_state.size() == S and index will not exceed bounds.
                return false; // Or throw an exception
            }
            int right=stack.back();
            stack.pop_back();
            int left=stack.back();
            stack.pop_back();
            
            bool parent_label = p_state[index++]; // Get the label for this internal node

            if(parent_label) // parent_label is true (1)
            {
                if(left==2 || right==2) return false; // '2' might represent 'any' or 'substituted'
                else stack.push_back(1);
            }
            else // parent_label is false (0)
            {
                if(left==0 && right==0) stack.push_back(0);
                else stack.push_back(2); // '2' representing a state that is not all zeros
            }
        }
    }
    return true;
}


// Check the 0-1 Constraint (if a node is 0, then at least one of its substree are all 0s)
// Parameter p is the StateKey representing the labeling of internal nodes.
bool Zero_One_Constraint(const StateKey& p_state) // p_state is the new name
{
    // S is the number of internal nodes, p_state.size() should be S.
    // StateKey is already in the correct order.

    vector<int> stack;
    int index = 0; // Index for accessing elements of p_state
    for(int i=0; i<trees[0].size(); i++) // Iterate through the tree structure string
    {
        if(trees[0][i]=='0' || trees[0][i]=='1')
        {
            int tmp=0;
            for(int j=0; j<N; j++) // N is the number of trees
                tmp+=trees[j][i]-'0';

            if(tmp>0) stack.push_back(1);
            else stack.push_back(0);
        }
        else // Encounter 'N', an internal node
        {
            if (index >= p_state.size()) {
                 // Error condition, similar to One_Oh_One_Constraint
                return false; 
            }
            int right=stack.back();
            stack.pop_back();
            int left=stack.back();
            stack.pop_back();

            bool parent_label_val = p_state[index++]; // Get the label (true for 1, false for 0)
            int current_parent_node_val_for_stack = parent_label_val ? 1 : 0;

            if(current_parent_node_val_for_stack == 0) // If parent label is 0
            {
                // If both subtrees are non-zero (i.e., >0), then it's a violation.
                // '0' means all-zero subtree. '1' means contains a '1'.
                // Original code: if(left>0 and right>0) return false;
                if(left > 0 && right > 0) return false; 
                // Original code: else if(left>0 or right>0) parent=1;
                // This means if one of them is >0 (contains a '1'), the effective value for stack is 1.
                // Otherwise (both are 0), it remains 0.
                else if(left > 0 || right > 0) current_parent_node_val_for_stack = 1;
            }
            // If parent_label_val was 1, current_parent_node_val_for_stack is already 1.
            stack.push_back(current_parent_node_val_for_stack);
        }
    }
    return true;
}

// Update the current tree
void UpdateCurrentTree(string tree, Tree* cur, Tree* pre)
{
    Valid_Internal_Labeling(tree); // This populates `validLabeling` with StateKey keys

    // S is the number of internal nodes, a member of TreeCentric.
    // This value (S) dictates the expected size of StateKey vectors.

    for(map<StateKey, int, StateKeyCompare>::iterator i = (pre->AccCost).begin(); i != (pre->AccCost).end(); ++i)
    {
        for(map<StateKey, int, StateKeyCompare>::iterator j = validLabeling.begin(); j != validLabeling.end(); ++j)
        {
            const StateKey& preV = i->first; // Accumulated labeling state from previous trees.
            int preCost = i->second;         // Cost associated with preV.
            const StateKey& curLabel_from_validLabeling = j->first; // Optimal labeling for the current tree.
            int curCost = j->second;         // Cost for curLabel_from_validLabeling.

            // Ensure StateKeys have the correct size (S, number of internal nodes).
            // This check is important as StateKey is now a dynamic vector.
            if (preV.size() != S || curLabel_from_validLabeling.size() != S) {
                cerr << "Error: StateKey size mismatch in UpdateCurrentTree." << endl;
                continue; 
            }

            // Calculate the new accumulated labeling state (curV) by ORing preV and curLabel.
            // This represents the combined state where a '1' appears if it was '1' in either previous accumulation or current tree's label.
            StateKey curV(S); 
            for(int k=0; k<S; ++k)
            {
                bool p_bit = preV[k];
                bool q_bit = curLabel_from_validLabeling[k];
                curV[k] = p_bit || q_bit; 
            }

            // If this combined state curV is new or offers a lower cost, update maps.
            if((cur->AccCost).count(curV)==0 || (cur->AccCost[curV]) > preCost + curCost)
            {
                (cur->AccCost)[curV] = preCost + curCost;
                // CurLabel in Tree stores the current tree's optimal label (curLabel_from_validLabeling) that leads to curV
                (cur->CurLabel)[curV] = curLabel_from_validLabeling; 
                // PreValue in Tree stores the previous accumulated state (preV) that leads to curV
                (cur->PreValue)[curV] = preV;
            }
        }
    }
}


TreeCentric::TreeCentric(vector<string> input) // Inline constructor definition
{
    trees=input;

    N=trees.size();
    if (N == 0 || trees[0].empty()) {
        S = 0;
        totalSubstitutions = 0;
        return;
    }
    S=trees[0].size()/2; // Number of internal nodes

    tree_layers.resize(N+1); // Use member variable tree_layers
    for(int i=0; i<N+1; i++) tree_layers[i]=new Tree(); // Use member variable tree_layers

    // Initialization
    if (S > 0) { 
        StateKey initialKey_AccCost(S, false); 
        (tree_layers[0]->AccCost)[initialKey_AccCost] = 0; // Use tree_layers
        StateKey dummyInitialCurLabel; 
        (tree_layers[0]->CurLabel)[initialKey_AccCost] = dummyInitialCurLabel; // Use tree_layers
    } else { 
        StateKey emptyKeyForS0;
        (tree_layers[0]->AccCost)[emptyKeyForS0] = 0; // Use tree_layers
        (tree_layers[0]->CurLabel)[emptyKeyForS0] = StateKey(); // Use tree_layers
    }

    for(int i=0; i<N; i++)
    {
        UpdateCurrentTree(trees[i], tree_layers[i+1], tree_layers[i]); // Use tree_layers
    }

    int totalSub = 1 << 30; 
    StateKey finalV;       
    bool finalV_found = false;

    if (S > 0) { 
        for(typename map<StateKey, int, StateKeyCompare>::iterator i = (tree_layers[N]->AccCost).begin(); i != (tree_layers[N]->AccCost).end(); ++i) // Use tree_layers
        {
            const StateKey& current_V = i->first;
            if (current_V.size() != S) {
                cerr << "Error: Key in AccCost has incorrect size before constraint checks." << endl;
                continue;
            }
            if(Zero_One_Constraint(current_V) && One_Oh_One_Constraint(current_V))
            {
                if(!finalV_found || (i->second) < totalSub)
                {
                    totalSub = i->second;
                    finalV = current_V;
                    finalV_found = true;
                }
            }
        }
    } else { 
        finalV = StateKey(); 
        totalSub = 0;
        finalV_found = true;
    }

    if (!finalV_found) {
        cerr << "Error: No optimal solution found that satisfies constraints." << endl;
        totalSubstitutions = -1; 
        optimalLabeling.assign(N, "Error: No solution");
        // The old delete loop for 'v' is removed. Destructor will handle tree_layers.
        return;
    }
    
    totalSubstitutions = totalSub;
    optimalLabeling.resize(N);

    if (S > 0) { 
        if ((tree_layers[N]->CurLabel).count(finalV)) { // Use tree_layers
             optimalLabeling[N-1] = printValue((tree_layers[N]->CurLabel).at(finalV)); // Use tree_layers
        } else {
            cerr << "Error: finalV not found in tree_layers[N]->CurLabel." << endl;
            optimalLabeling[N-1] = "Error: Label not found";
        }

        for(int i=N-1; i>0; i--)
        {
            if ((tree_layers[i+1]->PreValue).count(finalV)) { // Use tree_layers
                finalV = (tree_layers[i+1]->PreValue).at(finalV); // Use tree_layers
                if ((tree_layers[i]->CurLabel).count(finalV)) { // Use tree_layers
                    optimalLabeling[i-1] = printValue((tree_layers[i]->CurLabel).at(finalV)); // Use tree_layers
                } else {
                    cerr << "Error: finalV not found in tree_layers[i]->CurLabel during traceback." << endl;
                    optimalLabeling[i-1] = "Error: Label not found";
                    break; 
                }
            } else {
                 cerr << "Error: finalV not found in tree_layers[i+1]->PreValue during traceback." << endl;
                 for(int k=i-1; k>=0; --k) optimalLabeling[k] = "Error: Traceback failed";
                 break; 
            }
        }
    } else { 
        for(int i=0; i<N; ++i) {
            string stmp = "";
            for(int j=0; j<trees[i].size(); ++j)
                if(trees[i][j]!='N') stmp+=string(1,trees[i][j]);
            optimalLabeling[i] = stmp;
        }
    }

    if (S > 0) {
        for(int i=0; i<N; i++)
        {
            if (optimalLabeling[i].rfind("Error:", 0) == 0) continue;
            string internal_node_labels = optimalLabeling[i]; 
            int label_idx = 0;
            string stmp = "";
            for(int j=0; j<trees[i].size(); j++) {
                if(trees[i][j]!='N') {
                    stmp += string(1, trees[i][j]);
                } else {
                    if (label_idx < internal_node_labels.length()) {
                        stmp += string(1, internal_node_labels[label_idx++]);
                    } else {
                        cerr << "Error: Not enough labels for internal nodes in tree " << i << endl;
                        stmp += '?'; 
                    }
                }
            }
            optimalLabeling[i] = stmp;
        }
    }
    // The old delete loop for 'v' is removed. Destructor will handle tree_layers.
}

TreeCentric::~TreeCentric() // Inline destructor definition
{
    for(Tree* tree_ptr : tree_layers) {
        if (tree_ptr != nullptr) { 
            delete tree_ptr;
        }
    }
    tree_layers.clear(); 
}

};

