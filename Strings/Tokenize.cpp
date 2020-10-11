#include <string>
#include <vector>

using namespace std;

/**
 * Splits the input string into a vector of tokens
 * @param input The string to split
 * @param delimiters A string containing all the characters to consider as delimiters
 * @return A vector of tokens
 */
vector<string> tokenize(string input, string delimiters) {
	vector<string> tokens;
	string curr = "";
	for(int i = 0; i < input.size(); i++) {
		bool isDelimiter = false;
		for(int j = 0; j < delimiters.size(); j++) {
			if(input[i] == delimiters[j]) {
				isDelimiter = true;
				break;
			}
		}
		if(isDelimiter) {
			if(curr.size() > 0) {
				tokens.push_back(curr);
				curr = "";
			}
		} else {
			curr += input[i];
		}
	}
	if(curr.size() > 0) {
		tokens.push_back(curr);
		curr = "";
	}
	return tokens;
}