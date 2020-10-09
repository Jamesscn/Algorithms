#include <unordered_map>
#include <string>

using namespace std;

//NumberToBase("C8", 16, "0123456789ABCDEF") gives 200
unsigned long long int BaseToNumber(string input, unsigned long long int base, string alphabet) {
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