#include <unordered_map>

unordered_map<long long int, long long int> fibonacciMap({{0, 0}, {1, 1}, {2, 1}});

/**
 * Returns F(n) % mod in O(log(n))
 */
long long int fibonacci(long long int n, long long int mod) {
	if(fibonacciMap.count(n) == 0) {
		long long int k = n / 2;
		if(n % 2 == 0) {
			fibonacciMap[n] = (fibonacci(k, mod) * (2 * fibonacci(k + 1, mod) - fibonacci(k, mod) + mod)) % mod;
		} else {
			fibonacciMap[n] = (fibonacci(k, mod) * fibonacci(k, mod) + fibonacci(k + 1, mod) * fibonacci(k + 1, mod)) % mod;
		}
	}
	return fibonacciMap[n];
}