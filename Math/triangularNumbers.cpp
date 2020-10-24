
#include <math.h>
using namespace std;

/**
 * check if a number is triangular in constant time
 */
bool isTriangular(unsigned long long n){
    return ((float)sqrt(8*n+1) == floor(floor((float)sqrt(8*n+1))));
}