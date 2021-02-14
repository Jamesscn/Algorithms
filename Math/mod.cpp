/**
 * Calculates the positive and negative modulus of two numbers
 */
long long int mod(long long int a, long long int b) {
	return (a % b + b) % b;
}