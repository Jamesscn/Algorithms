/**
 * Calculates the greatest common divisor between two numbers
 */
long long gcd(long long a, long long b){
   return b == 0 ? a : gcd(b, a % b);
}