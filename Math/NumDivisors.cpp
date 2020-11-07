#include <vector>
#include <algorithm>

long long int NumDivisors(vector<long long int> primeFactors) {
	sort(primeFactors.begin(), primeFactors.end());
	long long int divisors = 1;
	long long int prevValue = 0;
	long long int count = 0;
	for(int i = 0; i < primeFactors.size(); i++) {
		if(primeFactors[i] == prevValue) {
			count++;
		} else {
			prevValue = primeFactors[i];
			divisors *= count + 1;
			count = 1;
		}
	}
	divisors *= count + 1;
	return divisors;
}