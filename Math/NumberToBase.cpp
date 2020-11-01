#include <string>

using namespace std;

/**
 * Converts an integer into a number in an arbitrary base stored in a string
 * For example, numberToBase(200, 16, "0123456789ABCDEF") returns "C8"
 * @param input An unsigned long integer to represent in the arbitrary base
 * @param base The arbitrary base
 * @param alphabet A list of characters representing each digit from zero to the base - 1
 * @return A string containing the representation of the number in the given base
 */
string numberToBase(unsigned long long int input, unsigned long long int base, string alphabet) {
	string output = "";
	while(input > 0) {
		output = alphabet[input % base] + output;
		input /= base;
	}
	return output;
}