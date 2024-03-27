#include <iostream>
#include <cmath>
#include <windows.h>

using namespace std;

const int N = pow(2,28);
int a[N];
int sum=0;
void init()
{
    for(int i =0;i<N;i++)
    {
        a[i]=i;
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
    for(int i=0;i<N;i++){
        sum+=a[i];
    }
     QueryPerformanceCounter(&tail);
    cout << "time: " << (tail.QuadPart - head.QuadPart) * 1000.0 / freq.QuadPart
         << "ms" << endl;
    return 0;
}
