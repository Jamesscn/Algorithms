#include <string>
#include <vector>

using namespace std;

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