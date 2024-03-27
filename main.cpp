#include <iostream>
#include <cmath>
#include <windows.h>

using namespace std;

const int N = pow(2,19);
int a[N];
int sum=0;
void init()
{
    for(int i =0;i<N;i++)
    {
        a[i]=i;
    }
}
int recursion(int n)
{
    if(n==1)return 0;

    else{
        for(int i=0;i<n/2;i++){
            a[i]+=a[n-i-1];
        }
        n=n/2;
        recursion(n);
    }
}
int main()
{
    LARGE_INTEGER head, tail, freq; // timers
    init();
     // similar to CLOCKS_PER_SEC
    QueryPerformanceFrequency(&freq);
    // start time
    QueryPerformanceCounter(&head);
    for(int i =0;i<100;i++){
    recursion(N);
    }
     QueryPerformanceCounter(&tail);
    cout << "time: " << (tail.QuadPart - head.QuadPart) * 1000.0 / freq.QuadPart
         << "ms" << endl;
    return 0;
}
