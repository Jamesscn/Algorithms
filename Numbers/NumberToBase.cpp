#include <string>

using namespace std;

//numberToBase(200, 16, "0123456789ABCDEF") gives "C8"
string numberToBase(unsigned long long int input, unsigned long long int base, string alphabet) {
	string output = "";
	while(input > 0) {
		output = alphabet[input % base] + output;
		input /= base;
	}
	return output;
}