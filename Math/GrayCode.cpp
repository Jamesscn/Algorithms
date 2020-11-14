/**
 * Converts a number N into gray code
 */
unsigned long long int GrayCode(unsigned long long int n) {
	return n ^ (n / 2);
}

/**
 * Converts a number G from gray code
 */
unsigned long long int InverseGrayCode(unsigned long long int g) {
	unsigned long long int n = 0;
	while(g != 0) {
		n ^= g;
		g /= 2;
	}
	return n;
}