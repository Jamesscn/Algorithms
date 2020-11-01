#include <vector>
#include <string>

using namespace std;

/**
 * Returns all permutations of str, including duplicates
 * By changing the vector into a set, the duplicates can be ommited
 * @param str The string to permute
 * @return A vector of permutations
 */
vector<string> Permute(string str) {
	vector<string> outputs;
	if(str.size() == 0) {
		outputs.push_back("");
		return outputs;
	}
	for(int i = 0; i < str.size(); i++) {
		string a = str.substr(0, i);
		string b = str.substr(i + 1, str.size() - i - 1);
		vector<string> subpermutations = Permute(a + b);
		for(int j = 0; j < subpermutations.size(); j++) {
			outputs.push_back(str[i] + subpermutations[j]);
		}
	}
	return outputs;
}