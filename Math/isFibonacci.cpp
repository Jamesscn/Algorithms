#include <math.h>

bool isSquare(unsigned long long x){
    unsigned long long s = sqrt(x);
    return (s*s == x);
}

bool isFib(unsigned long long n){
    return isSquare(5*n*n+4)||isSquare(5*n*n-4);
}