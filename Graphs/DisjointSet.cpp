#include <vector>

using namespace std;

/**
 * A class which allows the user to create a simple union-find/disjoint set data structure
 */
class DisjointSet {
	public:
	int numNodes;
	vector<int> parents;
	vector<vector<int>> sets;

	DisjointSet(int size) {
		numNodes = size;
		for(int i = 0; i < size; i++) {
			vector<int> currentSet = {i};
			sets.push_back(currentSet);
			parents.push_back(i);
		}
	}

	int find(int v) {
		return parents[v];
	}

	void join(int a, int b) {
		int firstSet = find(a);
		int secondSet = find(b);
		if(firstSet != secondSet) {
			if(sets[secondSet].size() < sets[firstSet].size()) {
				int tmp = firstSet;
				firstSet = secondSet;
				secondSet = tmp;
			}
			while(sets[firstSet].size() > 0) {
				int changeNode = sets[firstSet].back();
				parents[changeNode] = secondSet;
				sets[secondSet].push_back(changeNode);
				sets[firstSet].pop_back();
			}
		}
	}

	bool isSame(int a, int b) {
		return find(a) == find(b);
	}
};