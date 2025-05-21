#include <iostream>
#include <map>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <set>
#include <vector>     // For std::vector
#include <algorithm>  // For std::lexicographical_compare
#include <map>        // For std::map (though already implicitly used)
#include <iostream>   // For cout, cerr (already present)
// #include <limits> // For std::numeric_limits if MAXINT is replaced by it

#define MAXINT (1<<30) // Keep MAXINT for now, or replace with std::numeric_limits<int>::max()

using namespace std;

// Define Key Type and Comparator for NodeCentric algorithm.
// NodeStateKey represents the labeling state (0, 1, 2, or 3) for each of the N trees
// at a particular node in the species tree.
// Using std::vector<unsigned char> allows for N > 32 trees, overcoming the limit
// of a long long bitmask where each tree's state was packed into 2 bits.
// Each element in the vector corresponds to a tree's state.
using NodeStateKey = std::vector<unsigned char>; 

// Custom comparator for NodeStateKey to be used in std::map.
// Enables NodeStateKey instances to be used as keys in ordered maps.
struct NodeStateKeyCompare {
    bool operator()(const NodeStateKey& a, const NodeStateKey& b) const {
        // Lexicographical comparison is a standard way to compare sequences.
        return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
    }
};

class Node
{
	public:
		map<NodeStateKey, int, NodeStateKeyCompare> changes;
		map<NodeStateKey, NodeStateKey, NodeStateKeyCompare> leftV, rightV;
		Node* left;
		Node* right;
		Node(){left=NULL; right=NULL;}
		Node(const NodeStateKey& v_key) { changes[v_key]=0; left=right=NULL; }
};


class NodeCentric
{
	public:

//Store the results;
int totalSubstitutions;
vector<string> optimalLabeling;

//N: the number of trees
int N;
// long long checkZero; // Original
NodeStateKey checkZero_key; // New representation

//The input trees
vector<string> trees;

//Labeling Results
vector<NodeStateKey> results;

map<NodeStateKey, int, NodeStateKeyCompare> leftMap, rightMap;
map<NodeStateKey, NodeStateKey, NodeStateKeyCompare> leftMapNode, rightMapNode;

// map<long long, int> leftZero, rightZero; // Original, commented out
// map<long long, int> leftNonZeroMin, rightNonZeroMin; 
// map<long long, long long> leftNonZeroNode, rightNonZeroNode; 


// vp stores possible state transitions for parent and child nodes during the 'go' step.
// Indexed by tree number, it holds pairs of (parent_state, child_state).
// Max 50 trees assumed by original fixed-size array.
vector<pair<int,int> > vp[50]; 

// Converts a NodeStateKey (vector of states for N trees) into a string representation.
// N_val is the number of trees, expected to match n_key.size().
string printValue(const NodeStateKey& n_key, int N_val) 
{
    string s = "";
    // N_val should typically be equal to n_key.size().
    // The loop should iterate based on the actual size of n_key.
    for (int i = 0; i < n_key.size(); ++i) 
    {
        // Assuming n_key[i] directly stores the value 0-3.
        // The original prepended, so (N-1-i) for direct order, or iterate backwards.
        // Or, build and reverse, or prepend as original.
        s = string(1, '0' + n_key[i]) + s; 
    }
    return s;
}

// Helper function to check if a NodeStateKey (vector of states) consists of all zeros.
// Used to identify the "all-zero" state, analogous to a long long value of 0.
inline bool is_key_all_zeros(const NodeStateKey& key) {
    for (unsigned char val : key) {
        if (val != 0) return false;
    }
    return true;
}

// Helper function to check compatibility with the 'checkZero' constraint logic.
// The original 'checkZero' was a bitmask ...010101, meaning the LSB of each 2-bit state was 1.
// The condition `(key_long_long & checkZero_long_long) == 0` meant that for each tree,
// its state must be 0 (binary 00) or 2 (binary 10), i.e., the LSB of its state must be 0.
inline bool is_key_compatible_with_checkZero(const NodeStateKey& key) {
    for (unsigned char val : key) {
        // States 1 (01) or 3 (11) have their LSB as 1, thus are not compatible.
        if (val == 1 || val == 3) return false; 
    }
    return true;
}


// Checks if more than one tree has a state with LSB set to 1 (i.e., state is 1 or 3).
// This logic was adapted from the original `moreThanTwoOnes(long long p)`.
bool moreThanTwoOnes(const NodeStateKey& p_key) 
{
    int cnt = 0;
    for (unsigned char val : p_key) {
        // p&1 in original meant checking the LSB of the 2-bit state (e.g. 01 or 11).
        // So, if val is 1 or 3, its LSB is 1.
        if (val & 1) cnt++;
    }
    return cnt > 1;
}


void PostOrderTraversal(Node* cur, const NodeStateKey& value_key) // Changed parameter type
{
    //if(cur->left==NULL or cur->right==NULL) return; // Original comment
    
    if(cur==NULL) return;

    // Need to handle cases where value_key might not be in leftV or rightV if using [].
    // Using .find() and then .at() or direct access after check is safer.
    // Sticking to [] for now to mirror original, assuming keys will exist.
    auto it_left = cur->leftV.find(value_key);
    if (it_left != cur->leftV.end()) {
        PostOrderTraversal(cur->left, it_left->second);
    } else {
        // Handle error or assume a default if key not found, though original implies it would be.
        // For now, if not found, don't recurse on that path or use a default NodeStateKey().
        // This path indicates a potential issue if a key is expected but not found.
        // PostOrderTraversal(cur->left, NodeStateKey()); // Example: Pass empty key
    }

    auto it_right = cur->rightV.find(value_key);
    if (it_right != cur->rightV.end()) {
        PostOrderTraversal(cur->right, it_right->second);
    } else {
        // Similar handling for right child
        // PostOrderTraversal(cur->right, NodeStateKey()); // Example: Pass empty key
    }

    results.push_back(value_key); // results is now vector<NodeStateKey>
}

// int getBit(long long p, int i) // Old signature
int getBit(const NodeStateKey& p_key, int tree_idx) // New signature, i is tree_idx
{
    // p_key stores the state for each tree. p_key[tree_idx] is the state (0-3).
    if (tree_idx < p_key.size()) {
        return p_key[tree_idx];
    }
    // Handle error: tree_idx out of bounds
    // Return a default or throw an exception. For now, returning -1 (or similar error indicator).
    return -1; // Error indicator
}

// void go(Node* left, Node* right, Node* child) // Old signature
void go(Node* left_node, Node* right_node, Node* child_node) // New signature
{
    // N is a member of NodeCentric class
    if (left_node->changes.empty() || right_node->changes.empty() || child_node->changes.empty()) {
        // Handle error: one of the nodes has no states, cannot proceed.
        // This might indicate an issue upstream or a tree structure not expected.
        // For now, one could throw an exception, log an error, or simply return.
        // Depending on how vp is used later, not clearing/filling it might be problematic.
        // Clearing all vp[i] might be a safe default.
        for(int i=0; i<N; ++i) vp[i].clear();
        return;
    }

    // Get the first (and assumedly only, for leaves) NodeStateKey from each node's 'changes' map.
    // If these nodes are internal, 'changes' might have multiple entries, but original code
    // took .begin()->first. This implies it expects a specific key or any representative key.
    // This part needs to be robust if .changes can be empty or have many keys for non-leaf nodes.
    // For now, mirroring original logic:
    const NodeStateKey& lv_key = (left_node->changes).begin()->first;
    const NodeStateKey& rv_key = (right_node->changes).begin()->first;
    const NodeStateKey& cv_key = (child_node->changes).begin()->first;

    // Ensure keys are of expected size N.
    if (lv_key.size() != N || rv_key.size() != N || cv_key.size() != N) {
        // Handle error: key size mismatch.
        for(int i=0; i<N; ++i) vp[i].clear();
        return;
    }

    for(int i=0; i<N; i++) // Iterate through each tree
    {
        vp[i].clear();

        int lb = lv_key[i]; // State of current tree i in left child's key
        int rb = rv_key[i]; // State of current tree i in right child's key
        int cb = cv_key[i]; // State of current tree i in the 'child' param's key (which is one of left_node or right_node)

        if( lb==0 && rb==0 )
        {
            vp[i].push_back(make_pair(0,0));
        }
        else
        {
            if(cb==0)
            {
                vp[i].push_back(make_pair(1,0)); // Parent state 1, Child state 0
                vp[i].push_back(make_pair(2,0)); // Parent state 2, Child state 0
            }
            else // cb is 1 or 2 (or 3, if states can be 3)
            {
                vp[i].push_back(make_pair(1,1));
                vp[i].push_back(make_pair(2,1));
                vp[i].push_back(make_pair(2,2));
            }
        }
    }
}

map<NodeStateKey, int, NodeStateKeyCompare> leftZero, rightZero;
map<NodeStateKey, int, NodeStateKeyCompare> leftNonZeroMin, rightNonZeroMin;
map<NodeStateKey, NodeStateKey, NodeStateKeyCompare> leftNonZeroNode, rightNonZeroNode;


// Recursively generates all valid combinations of parent and child node states for the left child.
// 'current_tree_idx': The index of the current tree (0 to N-1) being processed.
// 'parentV_key': The NodeStateKey for the parent node being constructed in this recursive path.
// 'childV_key': The NodeStateKey for the child node (left child) being constructed.
// 'sub': Accumulated substitution cost for this path.
// 'child_node': Pointer to the actual left child Node object.
void allCombination(int current_tree_idx, NodeStateKey parentV_key, NodeStateKey childV_key, int sub, Node* child_node)
{
    // N is the total number of trees, a member of NodeCentric class.
    // parentV_key and childV_key are passed by value, ensuring each recursive path has its own state copy.
    if(current_tree_idx == N) // Base case: All N trees' states have been set for this combination.
    {
        // Check if the constructed childV_key is a valid state found in the child_node's 'changes' map.
        if((child_node->changes).count(childV_key))
        {
            // Apply 'checkZero' logic (now encapsulated in helper functions).
            // This logic partitions results based on parentV_key's compatibility and childV_key being all zeros.
            if(is_key_compatible_with_checkZero(parentV_key))
            {
                if(is_key_all_zeros(childV_key)) 
                    leftZero[parentV_key] = sub + (child_node->changes).at(childV_key); 
                else 
                {
                    if(leftNonZeroMin.count(parentV_key)==0 || leftNonZeroMin.at(parentV_key) > sub + (child_node->changes).at(childV_key))
                    {
                        leftNonZeroMin[parentV_key] = sub + (child_node->changes).at(childV_key);
                        leftNonZeroNode[parentV_key] = childV_key;
                    }
                }
            }
            else // parentV_key is not compatible with checkZero (i.e., has state 1 or 3 for some tree)
            {
                // Original condition: if( (childV & checkZero)==0 and childV!=0 ) return;
                // Translated: if child is compatible with checkZero AND child is not all-zeros, then this path is invalid.
                if( is_key_compatible_with_checkZero(childV_key) && !is_key_all_zeros(childV_key) ) return;

                // Otherwise, update the main leftMap for this parentV_key.
                if(leftMap.count(parentV_key)==0 || leftMap.at(parentV_key) > sub + (child_node->changes).at(childV_key))
                {
                    leftMap[parentV_key] = sub + (child_node->changes).at(childV_key);
                    leftMapNode[parentV_key] = childV_key;
                }
            }
        }
        return; // End of this recursive path.
    }
    else // Recursive step: Iterate through possible states for the current_tree_idx.
    {
        // vp[current_tree_idx] contains pairs of (parent_state_for_this_tree, child_state_for_this_tree).
        for(size_t i = 0; i < vp[current_tree_idx].size(); ++i)
        {
            unsigned char p_val = vp[current_tree_idx][i].first; 
            unsigned char q_val = vp[current_tree_idx][i].second;

            NodeStateKey next_parentV_key = parentV_key; 
            NodeStateKey next_childV_key = childV_key;
            
            // Ensure keys are sized correctly. This should ideally be guaranteed by initial call if N > 0.
            // If N=0, keys are empty and loops won't run, which is fine.
            if (next_parentV_key.size() != N && N > 0) next_parentV_key.resize(N); 
            if (next_childV_key.size() != N && N > 0) next_childV_key.resize(N);   

            if (N > 0) { // Only assign if N > 0, otherwise keys are empty.
              next_parentV_key[current_tree_idx] = p_val;
              next_childV_key[current_tree_idx] = q_val;
            }
            
            // Cost calculation: (p_val & 1) is LSB of parent's state, (q_val & 1) is LSB of child's state.
            // XOR means cost increases by 1 if LSBs are different.
            allCombination(current_tree_idx + 1, next_parentV_key, next_childV_key, sub + ((p_val & 1) ^ (q_val & 1)), child_node);
        }
    }
}

// Symmetric to allCombination, but for the right child and using rightMap, rightZero, etc.
void allCombination2(int current_tree_idx, NodeStateKey parentV_key, NodeStateKey childV_key, int sub, Node* child_node) 
{
    // Logic is symmetric to allCombination, just uses rightMap, rightZero, etc.
    if(current_tree_idx == N)
    {
        if((child_node->changes).count(childV_key))
        {
            if(is_key_compatible_with_checkZero(parentV_key))
            {
                if(is_key_all_zeros(childV_key))
                    rightZero[parentV_key] = sub + (child_node->changes).at(childV_key);
                else
                {
                    if(rightNonZeroMin.count(parentV_key)==0 || rightNonZeroMin.at(parentV_key) > sub + (child_node->changes).at(childV_key))
                    {
                        rightNonZeroMin[parentV_key] = sub + (child_node->changes).at(childV_key);
                        rightNonZeroNode[parentV_key] = childV_key;
                    }
                }
            }
            else
            {
                if( is_key_compatible_with_checkZero(childV_key) && !is_key_all_zeros(childV_key) ) return;

                if(rightMap.count(parentV_key)==0 || rightMap.at(parentV_key) > sub + (child_node->changes).at(childV_key))
                {
                    rightMap[parentV_key] = sub + (child_node->changes).at(childV_key);
                    rightMapNode[parentV_key] = childV_key;
                }
            }
        }
        return;
    }
    else
    {
        for(size_t i = 0; i < vp[current_tree_idx].size(); ++i)
        {
            unsigned char p_val = vp[current_tree_idx][i].first;
            unsigned char q_val = vp[current_tree_idx][i].second;

            NodeStateKey next_parentV_key = parentV_key;
            NodeStateKey next_childV_key = childV_key;
            if (next_parentV_key.size() != N) next_parentV_key.resize(N);
            if (next_childV_key.size() != N) next_childV_key.resize(N);

            next_parentV_key[current_tree_idx] = p_val;
            next_childV_key[current_tree_idx] = q_val;
            
            allCombination2(current_tree_idx + 1, next_parentV_key, next_childV_key, sub + ((p_val & 1) ^ (q_val & 1)), child_node);
        }
    }
}

NodeCentric(vector<string> input)
{
    trees = input;
    N = trees.size();
    if (N == 0) { // Handle empty input
        totalSubstitutions = 0;
        return;
    }
    // Initialize checkZero_key. Its direct value isn't used by helper functions,
    // but it's initialized to reflect the change from long long.
    // The core logic is now within is_key_compatible_with_checkZero.
    // Original checkZero was a bitmask ...010101.
    checkZero_key.assign(N, 1); // Assigns N elements, each with value 1.

    vector<Node*> stack; // Stack to build the species tree representation.

    // Iterate through the postfix representation of the species tree.
    for (size_t i = 0; i < trees[0].size(); ++i) 
    {
        if (trees[0][i] == '0' || trees[0][i] == '1') // Leaf node in species tree string
        {
            // Construct the initial NodeStateKey for this leaf.
            // This key represents the observed states of this leaf across all N gene trees.
            NodeStateKey current_leaf_state(N);
            for (int j = 0; j < N; ++j) // For each gene tree
            {
                if (i < trees[j].size()) { 
                    current_leaf_state[j] = trees[j][i] - '0'; // State is '0' or '1'
                } else {
                    // Should not happen with valid input where all tree strings are consistent.
                    current_leaf_state[j] = 0; // Default/error value
                }
            }
            Node* newNode = new Node(current_leaf_state); // Create a new tree Node.
            stack.push_back(newNode);
        }
        else // Internal node ('N') in species tree string
        {
            if (stack.size() < 2) {
                // Error: Not enough operands for internal node
                // Handle this error appropriately, e.g., by cleaning up and returning.
                cerr << "Error: Tree structure implies internal node but stack has < 2 elements." << endl;
                for(Node* n_ptr : stack) delete n_ptr; // Clean up allocated nodes
                totalSubstitutions = -1; // Indicate error
                return;
            }
            Node* right_node = stack.back(); stack.pop_back();
            Node* left_node = stack.back(); stack.pop_back();

            Node* newNode = new Node(); // Default constructor
            newNode->left = left_node;
            newNode->right = right_node;

            leftMap.clear(); rightMap.clear();
            leftMapNode.clear(); rightMapNode.clear();
            leftZero.clear(); rightZero.clear();
            leftNonZeroMin.clear(); rightNonZeroMin.clear();
            leftNonZeroNode.clear(); rightNonZeroNode.clear();

            go(left_node, right_node, left_node); // vp for left child combinations
            allCombination(0, NodeStateKey(N,0), NodeStateKey(N,0), 0, left_node); // Initial keys are all zeros

            go(left_node, right_node, right_node); // vp for right child combinations
            allCombination2(0, NodeStateKey(N,0), NodeStateKey(N,0), 0, right_node); // Initial keys are all zeros

            // Combine results from left and right children
            // Iterate through non-zero combinations (leftMap, rightMap)
            for (map<NodeStateKey, int, NodeStateKeyCompare>::iterator it = leftMap.begin(); it != leftMap.end(); ++it)
            {
                const NodeStateKey& value_key = it->first;
                if (rightMap.count(value_key)) {
                    (newNode->changes)[value_key] = it->second + rightMap.at(value_key);
                    (newNode->leftV)[value_key] = leftMapNode.at(value_key);
                    (newNode->rightV)[value_key] = rightMapNode.at(value_key);
                }
            }

            NodeStateKey zero_state_key(N, 0); // Represents the "all zeros" state key

            // Handle combinations involving zero states
            for (map<NodeStateKey, int, NodeStateKeyCompare>::iterator it = leftZero.begin(); it != leftZero.end(); ++it)
            {
                const NodeStateKey& value_key = it->first;
                if (rightZero.count(value_key)) {
                    int current_cost = it->second + rightZero.at(value_key);
                    if ((newNode->changes).count(value_key) == 0 || (newNode->changes).at(value_key) > current_cost) {
                        (newNode->changes)[value_key] = current_cost;
                        (newNode->leftV)[value_key] = zero_state_key;
                        (newNode->rightV)[value_key] = zero_state_key;
                    }
                }
                if (rightNonZeroMin.count(value_key)) {
                    int current_cost = it->second + rightNonZeroMin.at(value_key);
                    if ((newNode->changes).count(value_key) == 0 || (newNode->changes).at(value_key) > current_cost) {
                        (newNode->changes)[value_key] = current_cost;
                        (newNode->leftV)[value_key] = zero_state_key;
                        (newNode->rightV)[value_key] = rightNonZeroNode.at(value_key);
                    }
                }
            }
    
            for (map<NodeStateKey, int, NodeStateKeyCompare>::iterator it = rightZero.begin(); it != rightZero.end(); ++it)
            {
                const NodeStateKey& value_key = it->first;
                // Skip if already handled by leftZero x rightZero combination
                // if (leftZero.count(value_key)) { /* already handled or would be symmetric */ }
                if (leftNonZeroMin.count(value_key)) {
                    int current_cost = leftNonZeroMin.at(value_key) + it->second;
                    if ((newNode->changes).count(value_key) == 0 || (newNode->changes).at(value_key) > current_cost) {
                        (newNode->changes)[value_key] = current_cost;
                        (newNode->leftV)[value_key] = leftNonZeroNode.at(value_key);
                        (newNode->rightV)[value_key] = zero_state_key;
                    }
                }
            }
            stack.push_back(newNode);
        }
    }

    if (stack.empty()) {
         cerr << "Error: Stack is empty at the end of processing. No root node." << endl;
         totalSubstitutions = -1;
         // Clean up any nodes if necessary (though stack is empty, other allocations might exist if error handling is more complex)
         return;
    }
    if (stack.size() > 1) {
        cerr << "Error: Stack has more than one node at the end. Tree structure malformed or error in processing." << endl;
        // Clean up all nodes on stack
        for(Node* n_ptr : stack) delete n_ptr;
        totalSubstitutions = -1;
        return;
    }

    // Find the optimal value from the root node's changes
    int minSub = MAXINT;
    NodeStateKey optimal_label_key; 
    bool found_optimal = false;

    if (!(stack[0]->changes).empty()) {
       for (map<NodeStateKey, int, NodeStateKeyCompare>::iterator it = (stack[0]->changes).begin(); it != (stack[0]->changes).end(); ++it)
       {
           if (it->second < minSub)
           {
               minSub = it->second;
               optimal_label_key = it->first;
               found_optimal = true;
           }
       }
    }
    
    if (!found_optimal) {
        // This case can happen if root node's 'changes' map is empty, or no valid label was found.
        // Set to a default state or error. Original code used label=-1 (long long).
        // For NodeStateKey, an empty key or a key of N zeros might be a default.
        // If minSub is still MAXINT, it means no solution.
        if (minSub == MAXINT) {
            cerr << "Warning: No optimal solution found. totalSubstitutions will be MAXINT." << endl;
            totalSubstitutions = MAXINT; // Or some other error indicator like -1
            // optimalLabeling will be empty or indicate error.
        }
        // If N > 0 and optimal_label_key is still default-constructed (empty), resize it.
        if (N > 0 && optimal_label_key.empty()) {
            optimal_label_key.assign(N, 0); // Default to all zeros if no label found but N > 0.
        }
    }
    totalSubstitutions = minSub;

    // Post-order Traversal to reconstruct results for each node
    results.clear(); // Clear previous results if any
    if (found_optimal) { // Only traverse if a valid optimal label was found
        PostOrderTraversal(stack[0], optimal_label_key);
    } else if (N > 0 && !stack.empty()) { // Try with a default key if nothing optimal found but tree exists
        PostOrderTraversal(stack[0], NodeStateKey(N,0)); // e.g. traverse with all-zero key
    }


    // Print labeling of each tree
    optimalLabeling.resize(N); // Ensure optimalLabeling is sized for N trees
    for (int i = 0; i < N; ++i) // Iterate for each tree string to construct
    {
        stringstream sss;
        for (size_t j = 0; j < results.size(); ++j) // Iterate through nodes in post-order
        {
            // results[j] is NodeStateKey for node j. results[j][i] is state of tree i at node j.
            // Original: ((results[j]>>(2*i))&1) - this extracts the LSB of the 2-bit state for tree i.
            // Example: if state is 0 (00), LSB is 0. If 1 (01), LSB is 1. If 2 (10), LSB is 0. If 3 (11), LSB is 1.
            if (i < results[j].size()) { // Check bounds for tree index i
                 sss << (results[j][i] & 1);
            } else {
                 // Should not happen if all NodeStateKeys in results have size N
                 sss << '0'; // Default or error character
            }
        }
        optimalLabeling[i] = sss.str(); // Store the string for tree i
    }
    
    // Final cleanup of the root node (and thus the whole tree)
    if (!stack.empty()) {
        delete stack[0];
        stack.clear();
    }
}

};
