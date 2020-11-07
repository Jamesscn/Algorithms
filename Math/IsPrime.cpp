/**
 * Checks if a value is prime in O(log^3 N) time using a deterministic implementation of the Miller-Rabin algorithm
 * This code only works for values less than 2^64, and requires __int128
 * Written by Brett Hale (StackOverflow 24096332)
 * @param n The value that will be tested
 * @return True if n is prime, and false if n is composite
 */
bool IsPrime (unsigned long long int n) {
	if(n == 1) {
		return 0;
	}
	const unsigned long long int sprp32_base[] = {2, 7, 61, 0};
	const unsigned long long int sprp64_base[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022, 0};
	const unsigned long long int *sprp_base;
	if((n & 1) == 0) {
		return n == 2;
	}
	sprp_base = (n <= 4294967295U) ? sprp32_base : sprp64_base;
	for(; *sprp_base != 0; sprp_base++) {
		unsigned long long int a = *sprp_base, m = n - 1, r, y, s = 1;
		while ((m & (1UL << s)) == 0) {
			s++;
		}
		r = m >> s;
		if ((a %= n) == 0) {
			continue;
		}
		unsigned __int128 u = 1, w = a;
		while (r != 0) {
			if ((r & 1) != 0) {
				u = (u * w) % n;
			}
			if ((r >>= 1) != 0) {
				w = (w * w) % n;
			}
		}
		if ((y = (unsigned long long int)u) == 1) {
			continue;
		}
		for (unsigned long long int j = 1; j < s && y != m; j++) {
			u = y;
			u = (u * u) % n;
			if ((y = (unsigned long long int)u) <= 1) {
				return false;
			}
		}
		if(y != m) {
			return false;
		}
	}
	return true;
}