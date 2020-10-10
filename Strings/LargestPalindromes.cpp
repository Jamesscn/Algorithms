#include <unordered_set>
#include <utility>
#include <string>

using namespace std;

/**
 * Returns all the largest palindromes in the string
 * @param input The input string
 * @return A pair which contains an int (the length of the largest palindrome) and an unordered set of these palindromes
 */
pair<int, unordered_set<string>> getLargestPalindromes(string input) {
	if(input.size() < 2) {
		unordered_set<string> palindrome;
		palindrome.insert(input);
		return make_pair(input.size(), palindrome);
	}
	if(input[0] == input[input.size() - 1]) {
		pair<int, unordered_set<string>> palindromeInfo = getLargestPalindromes(input.substr(1, input.size() - 2));
		palindromeInfo.first += 2;
		unordered_set<string> newSet;
		for(const string &str: palindromeInfo.second) {
			newSet.insert(input[0] + str + input[input.size() - 1]);
		}
		palindromeInfo.second = newSet;
		return palindromeInfo;
	}
	pair<int, unordered_set<string>> palindromeInfo;
	pair<int, unordered_set<string>> a = getLargestPalindromes(input.substr(0, input.size() - 1));
	pair<int, unordered_set<string>> b = getLargestPalindromes(input.substr(1, input.size() - 1));
	palindromeInfo.first = max(a.first, b.first);
	if(a.first == palindromeInfo.first) {
		for(const string &str: a.second) {
			palindromeInfo.second.insert(str);
		}
	}
	if(b.first == palindromeInfo.first) {
		for(const string &str: b.second) {
			palindromeInfo.second.insert(str);
		}
	}
	return palindromeInfo;
}