void fastscan(int &x){
   char ch; bool f= 0; int a=0;
   while(!((((ch=getchar())>='0')&&(ch<='9'))||(ch=='-')));
   if(ch!='-')a*=10,a+=ch-'0';else f=1;
   while(((ch=getchar())>='0')&&(ch<='9'))a*=10, a+=ch-'0';
   if(f)a=-a;x=a;
}