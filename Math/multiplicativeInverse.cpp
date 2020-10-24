#define mod 1000000007

long long inve(long long a){
   long long b = mod-2, ans = 1;
   while (b) {
       if (b&1) {
           ans = (ans * a) % mod;
       }
       a = (a * a) % mod;
       b >>= 1;
   }
   return ans;
}