#include <vector>

long long int Brent(long long int N) {
	if(N % 2 == 0) {
		return 2;
	}
	long long int y = (rand() % (N - 1)) + 1;
	long long int c = (rand() % (N - 1)) + 1;
	long long int m = (rand() % (N - 1)) + 1;
	long long int x, k, ys, g = 1, r = 1, q = 1;
	while(g == 1) {            
		x = y;
		for(int i = 0; i < r; i++) {
			y = (modMult(y, y, N) + c) % N;
		}
		k = 0;
		while(k < r and g == 1) {
			ys = y;
			for(int i = 0; i < min(m, r - k); i++) {
				y = (modMult(y, y, N) + c) % N;
				q = modMult(q, abs(x - y), N);
			}
			g = gcd(q, N);
			k = k + m;
		}
		r = r * 2;
	}
	if(g == N) {
		while(true) {
			ys = (modMult(ys, ys, N) + c) % N;
			g = gcd(abs(x - ys), N);
			if(g > 1) {
				break;
			}
		}
	}
	return g;
}

/**
 * Finds all prime factors of the given value using Brent's method
 * A fast IsPrime algorithm like Miller-Rabin's deterministic algorithm should be used
 * Disclamer: factors are not always in order
 */
vector<long long int> PrimeFactorization(long long int value) {
	vector<long long int> factors;
	if(value < 2) {
		return factors;
	}
	vector<long long int> factorsToCheck;
	factorsToCheck.push_back(value);
	while(factorsToCheck.size() > 0) {
		long long int currentFactor = factorsToCheck.back();
		factorsToCheck.pop_back();
		if(currentFactor == 1) {
			continue;
		}
		if(IsPrime(currentFactor)) {
			factors.push_back(currentFactor);
		} else {
			long long int newFactor = Brent(currentFactor);
			factorsToCheck.push_back(newFactor);
			factorsToCheck.push_back(currentFactor / newFactor);
		}
	}
	return factors;
}