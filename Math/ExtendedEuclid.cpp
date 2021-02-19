#include <utility>

using namespace std;

/**
 * Applies the extended euclid algorithm to find a solution to the equation a*x + b*y = gcd(a, b).
 * @param a The first value.
 * @param b The second value.
 * @return A pair of long long ints with the solutions x and y found by the algorithm.
 */
pair<long long int, long long int> ExtendedEuclid(long long int a, long long int b) {
    if (b == 0) {
        return make_pair(1, 0);
    }
    pair<long long int, long long int> result = ExtendedEuclid(b, a % b);
    return make_pair(result.second, result.first - result.second * (a / b));
}