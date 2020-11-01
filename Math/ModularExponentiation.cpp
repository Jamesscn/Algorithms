/**
 * Calculates the operation base ^ exponent % modulus for long integers
 */
long long int modularExponentiation(long long int base, long long int exponent, long long int modulus) {
	long long int power = 1;
	base %= modulus;
	while(exponent > 0) {
		if(exponent % 2) {
			power = (power * base) % modulus;
		}
		exponent /= 2;
		base = (base * base) % modulus;
	}
	return power;
}