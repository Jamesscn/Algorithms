#include <unordered_map>
#include <string>

using namespace std;

/**
 * Converts a number in an arbitrary base stored in a string into an actual integer
 * For example, baseToNumber("C8", 16, "0123456789ABCDEF") returns 200
 * @param input A string containing a number in an arbitrary base
 * @param base The arbitrary base
 * @param alphabet A list of characters representing each digit from zero to the base - 1
 * @return An unsigned long integer which is equivalent to the input number
 */
unsigned long long int baseToNumber(string input, unsigned long long int base, string alphabet) {
	unsigned long long int output = 0;
	unordered_map<char, unsigned long long int> alphabetMap;
	for(unsigned long long int i = 0; i < alphabet.size(); i++) {
		alphabetMap[alphabet[i]] = i;
	}
	for(unsigned long long int i = 0; i < input.size(); i++) {
		output *= base;
		output += alphabetMap[input[i]];
	}
	return output;
}