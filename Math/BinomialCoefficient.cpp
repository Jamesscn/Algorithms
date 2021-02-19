#include <algorithm>
#include <cmath>

using namespace std;

/**
 * Finds the binomial coefficient for a given n and r in O(r).
 * @param n The top value of the binomial coefficient.
 * @param r The bottom value of the binomial coefficient.
 */
long long int BinomialCoefficient(int n, int r) {
    long double coefficient = 1;
    for(int i = 1; i <= r; i++) {
        coefficient *= (long double)(n - r + i) / i;
	}
    return round(coefficient);
}