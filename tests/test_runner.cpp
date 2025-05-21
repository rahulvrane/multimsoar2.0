#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <utility> // For std::pair
#include <tuple>   // For std::tuple
#include <cassert>
#include <limits>  // For std::numeric_limits

// It might be necessary to include the .cpp files directly if headers alone
// are not sufficient due to template instantiations or if not using a build system
// that compiles and links them separately. This is common for simple test runners.
// However, try with headers first.
#include "../src/Hungarian.h"
// For Partition_modified, we need its definition and other related types/functions
// This might require including MultiMSOARSoftware.cpp or refactoring parts of it
// for testability. For this subtask, let's focus on Hungarian.h which is self-contained.
// If testing Partition_modified is too complex due to dependencies from MultiMSOARSoftware.cpp
// in this environment, stick to Hungarian.h tests.

// --- Test for Hungarian Algorithm ---
void test_hungarian_basic() {
    std::cout << "Running test_hungarian_basic..." << std::endl;
    Hungarian::matrix weights1 = {{10, 5, 2}, {7, 8, 3}, {6, 0, 9}};
    Hungarian h1(weights1);
    // Expected matching: (0,0)=10, (1,1)=8, (2,2)=9. Total = 27
    // matchingX[0] should be 0, matchingX[1] should be 1, matchingX[2] should be 2.
    assert(h1.totalweight == 27);
    assert(h1.matchingX[0] == 0);
    assert(h1.matchingX[1] == 1);
    assert(h1.matchingX[2] == 2);
    std::cout << "test_hungarian_basic PASSED" << std::endl;
}

void test_hungarian_different_optimal() {
    std::cout << "Running test_hungarian_different_optimal..." << std::endl;
    Hungarian::matrix weights2 = {{1, 2, 8}, {3, 7, 4}, {6, 5, 0}};
    Hungarian h2(weights2);
    // Optimal: (0,2)=8, (1,1)=7, (2,0)=6. Total = 21.
    // matchingX[0]=2, matchingX[1]=1, matchingX[2]=0
    assert(h2.totalweight == 21);
    assert(h2.matchingX[0] == 2);
    assert(h2.matchingX[1] == 1);
    assert(h2.matchingX[2] == 0);
    std::cout << "test_hungarian_different_optimal PASSED" << std::endl;
}

void test_hungarian_one_by_one() {
    std::cout << "Running test_hungarian_one_by_one..." << std::endl;
    Hungarian::matrix weights = {{100}};
    Hungarian h(weights);
    assert(h.totalweight == 100);
    assert(h.matchingX[0] == 0);
    std::cout << "test_hungarian_one_by_one PASSED" << std::endl;
}

// --- Placeholder for future tests ---
// void test_partition_logic() { ... }
// void test_treelabeling_logic() { ... }


int main(int argc, char** argv) {
    std::cout << "--- Starting MultiMSOAR Unit Tests ---" << std::endl;
    
    test_hungarian_basic();
    test_hungarian_different_optimal();
    test_hungarian_one_by_one();

    // Add calls to other test functions here
    // e.g., test_partition_logic();

    std::cout << "--- All MultiMSOAR Unit Tests Completed ---" << std::endl;
    return 0;
}
