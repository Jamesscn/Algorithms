/**
 * Calculates and stores all primes less than or equal to the limit in O(N) time
 * @param primes An array where all the primes less than or equal to the limit shall be stored
 * @param limit The largest value that will be checked by the function
 * @return The amount of primes that are stored in the array of primes
 */
int SieveOfAtkin(long long int* primes, long long int limit) {
    int index = 0;
    limit++;
    if(limit < 3) {
        return index;
    }
    primes[index] = 2;
    index++;
    if(limit < 4) {
        return index;
    }
    primes[index] = 3;
    index++;
    bool* sieve = new bool[limit];
    for(long long int i = 0; i < limit; i++) {
        sieve[i] = false;
    }
    for(long long int i = 1; i * i < limit; i++) {
        for(long long int j = 1; j * j < limit; j++) {
            long long int a = 4 * i * i + j * j;
            long long int b = 3 * i * i + j * j;
            long long int c = 3 * i * i - j * j;
            if(a <= limit && (a % 12 == 1 || a % 12 == 5)) {
                sieve[a] = !sieve[a];
            }
            if(b <= limit && b % 12 == 7) {
                sieve[b] = !sieve[b];
            }
            if(c <= limit && i > j && c % 12 == 11) {
                sieve[c] = !sieve[c];
            }
        }
    }
    for(long long int i = 5; i * i < limit; i++) {
        if(sieve[i]) {
            for(long long int j = i * i; j < limit; j += i * i) {
                sieve[j] = false;
            }
        }
    }
    for(long long int i = 5; i < limit; i++) {
        if(sieve[i]) {
            primes[index] = i;
            index++;
        }
    }
    return index;
}