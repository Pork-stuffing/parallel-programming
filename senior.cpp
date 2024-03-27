#include <iostream>
#include <windows.h>
#include <stdlib.h>

using namespace std;

const int N = 4000; // matrix size
double b[N][N], sum[N],a[N];

void init(int n) // generate a N*N matrix
{
    for (int i = 0; i < N; i++)
    {
        a[i]=2*i;
        for (int j = 0; j < N; j++)
            b[i][j] = i + j;
    }
}
//cacheÓÅ»¯Ëã·¨
int main()
{
    LARGE_INTEGER head, tail, freq; // timers
    init(N);
    // similar to CLOCKS_PER_SEC
    QueryPerformanceFrequency(&freq);
    // start time
    QueryPerformanceCounter(&head);
    for (int j = 0; j < N; j++)
    {
        for (int i = 0; i < N; i++)
            sum[i] += b[j][i]*a[i];
    }

    // end time
    QueryPerformanceCounter(&tail);
    cout << "time: " << (tail.QuadPart - head.QuadPart) * 1000.0 / freq.QuadPart
         << "ms" << endl;
    return 0;
}
