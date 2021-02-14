/**
 * Calculates the least common multiple between two numbers
 */
long long gcd(long long a, long long b){
   return b == 0 ? a : gcd(b, a % b);
}

long long lcm(long long a, long long b){
   return (a * b) / gcd(a, b);
}