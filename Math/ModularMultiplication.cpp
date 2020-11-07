/**
 * Calculates the operation a * b % modulus for unsigned long integers without overflow
 */
unsigned long long int modMult(unsigned long long int a, unsigned long long int b, unsigned long long int modulus) {
	unsigned long long int product = 0, tempB;
	if(b >= modulus) {
		if(modulus > 18446744073709551615ULL / 2ULL) {
			b -= modulus;
		} else {
			b %= modulus;
		}
	}
	while(a != 0) {
		if(a % 2 == 1) {
			if(b >= modulus - product) {
				product -= modulus;
			}
			product += b;
		}
		a /= 2;
		tempB = b;
		if(b >= modulus - b) {
			tempB -= modulus;
		}
		b += tempB;
	}
	return product;
}