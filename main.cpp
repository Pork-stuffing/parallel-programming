#include <iostream>
#include<emmintrin.h>
#include<time.h>
#include<Windows.h>
#include<math.h>
#include <immintrin.h>
using namespace std;

const int N = 2000;
float elm[N][N]={0};

void reset(float **test)
{
    for(int i=0;i<N;i++)
    {
        for(int j = 0;j<N;j++)
        {
            test[i][j]=elm[i][j];
        }
    }
    return;
}

//ƽ���㷨
void simple(float **e)
{
    for(int k=0;k<N;k++)
    {
        for(int j=k;j<N;j++)
        {
            e[k][j]/=e[k][k];
        }
        for(int i = k+1;i<N;i++)
        {
            for(int j = k+1;j<N;j++)
            {
                e[i][j] = e[i][j] - e[i][k]*e[k][j];
            }
            e[i][k]=0;
        }
    }
    return;
}

//SSE�㷨
void sse_gausseliminate(float** A)//����
{
    __m128 t1, t2, t3;//��λ�����ȹ��ɵ�����
    for (int k = 0; k < N; k++)
    {
        int preprocessnumber = (N - k - 1) % 4;//Ԥ���������,�ܱ�������
        int begin = k + 1 + preprocessnumber;
        __attribute__((aligned(16)))float head[4] = { A[k][k],A[k][k],A[k][k],A[k][k] };
        t2 = _mm_load_ps(head);
        for (int j = k + 1; j < k + 1 + preprocessnumber; j++)
        {
            A[k][j] = A[k][j] / A[k][k];
        }
        for (int j = begin; j < N; j += 4)
        {
            t1 = _mm_load_ps(A[k] + j);
            t1 = _mm_div_ps(t1, t2);
            _mm_store_ss(A[k] + j, t1);
        }
        A[k][k] = 1;
        t1 = _mm_setzero_ps();//����
        t2 = _mm_setzero_ps();
        //��ȥͷ��Ϊ���ĸ��ĸ��Ĵ���
        for (int i = k + 1; i < N; i++)
        {
            for (int j = k + 1; j < k + 1 + preprocessnumber; j++)
            {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
            A[i][k] = 0;
        }
        for (int i = k + 1; i < N; i++)
        {
            __attribute__((aligned(16)))float head1[4] = { A[i][k],A[i][k],A[i][k],A[i][k] };
            t3 = _mm_load_ps(head1);
            for (int j = begin; j < N; j += 4)
            {
                t1 = _mm_load_ps(A[k] + j);
                t2 = _mm_load_ps(A[i] + j);
                t1 = _mm_mul_ps(t1, t3);
                t2 = _mm_sub_ps(t2, t1);
                _mm_store_ss(A[i] + j, t2);
            }
            A[i][k] = 0;
        }
    }
}

//�������SSE�㷨
void SIMD_notaligned_SSE_gausseliminate(float** A)
{
    __m128 t1, t2, t3;//��λ�����ȹ��ɵ�����
    for (int k = 0; k < N; k++)
    {
        int preprocessnumber = (N - k - 1) % 4;//Ԥ���������,�ܱ�������
        int begin = k + 1 + preprocessnumber;
        float head[4] = { A[k][k],A[k][k],A[k][k],A[k][k] };
        t2 = _mm_loadu_ps(head);
        for (int j = k + 1; j < k + 1 + preprocessnumber; j++)
        {
            A[k][j] = A[k][j] / A[k][k];
        }
        for (int j = begin; j < N; j += 4)
        {
            t1 = _mm_loadu_ps(A[k] + j);
            t1 = _mm_div_ps(t1, t2);
            _mm_store_ss(A[k] + j, t1);
        }
        A[k][k] = 1;
        t1 = _mm_setzero_ps();//����
        t2 = _mm_setzero_ps();
        //��ȥͷ��Ϊ���ĸ��ĸ��Ĵ���
        for (int i = k + 1; i < N; i++)
        {
            for (int j = k + 1; j < k + 1 + preprocessnumber; j++)
            {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
            A[i][k] = 0;
        }
        for (int i = k + 1; i < N; i++)
        {
            float head1[4] = { A[i][k],A[i][k],A[i][k],A[i][k] };
            t3 = _mm_loadu_ps(head1);
            for (int j = begin; j < N; j += 4)
            {
                t1 = _mm_loadu_ps(A[k] + j);
                t2 = _mm_loadu_ps(A[i] + j);
                t1 = _mm_mul_ps(t1, t3);
                t2 = _mm_sub_ps(t2, t1);
                _mm_store_ss(A[i] + j, t2);
            }
            A[i][k] = 0;
        }
    }
}

//AVX�㷨
void avx_gauss(float** e)
{
    __m256 t1, t2, t3;
    for (int k = 0; k < N; k++)
    {
        float temp1[8] = { e[k][k], e[k][k], e[k][k], e[k][k], e[k][k], e[k][k], e[k][k], e[k][k] };
        t1 = _mm256_loadu_ps(temp1);
        int j = k + 1;
        for (j; j < N - 7; j += 8)
        {
            t2 = _mm256_loadu_ps(e[k] + j);
            t3 = _mm256_div_ps(t2, t1);
            _mm256_storeu_ps(e[k] + j, t3);
        }
        for (j; j < N; j++)
            e[k][j] = e[k][j] / e[k][k];

        e[k][k] = 1.0;

        for (int i = k + 1; i < N; i++)
        {
            float temp2[8] = { e[i][k], e[i][k], e[i][k], e[i][k], e[i][k], e[i][k], e[i][k], e[i][k] };
            t1 = _mm256_loadu_ps(temp2);
            j = k + 1;
            for (j; j < N - 7; j += 8)
            {
                t2 = _mm256_loadu_ps(e[i] + j);
                t3 = _mm256_loadu_ps(e[k] + j);
                t3 = _mm256_mul_ps(t1, t3);
                t2 = _mm256_sub_ps(t2, t3);
                _mm256_storeu_ps(e[i] + j, t2);
            }
            for (j; j < N; j++)
                e[i][j] = e[i][j] - e[i][k] * e[k][j];

            e[i][k] = 0;
        }
    }
}

int main()
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            elm[i][j] = (rand() % 100);
        }
    }
    float** test = new float*[N];
    for (int i = 0; i < N; i++)
    {
        test[i] = new float[N];
    }

    reset(test);

     srand(time(NULL));
    LARGE_INTEGER timeStart;	//��ʼʱ��
    LARGE_INTEGER timeEnd;		//����ʱ��

    LARGE_INTEGER frequency;	//��ʱ��Ƶ��
    QueryPerformanceFrequency(&frequency);
    double quadpart = (double)frequency.QuadPart;//��ʱ��Ƶ��

    //ƽ���㷨
    QueryPerformanceCounter(&timeStart);
    simple(test);
    QueryPerformanceCounter(&timeEnd);
    double _Simple = (timeEnd.QuadPart - timeStart.QuadPart) / quadpart;
    printf("Simple:%f\n", _Simple);
    cout << endl;
    reset(test);

    //�����SSE�㷨
    QueryPerformanceCounter(&timeStart);
    sse_gausseliminate(test);
    QueryPerformanceCounter(&timeEnd);
    double _SSE_Gauss = (timeEnd.QuadPart - timeStart.QuadPart) / quadpart;
    printf("SSE_Gauss:%f\n", _SSE_Gauss);
    cout << endl;
    reset(test);

    //�������SSE�㷨
    QueryPerformanceCounter(&timeStart);
    SIMD_notaligned_SSE_gausseliminate(test);
    QueryPerformanceCounter(&timeEnd);
    double notaligned_SSE_Gauss = (timeEnd.QuadPart - timeStart.QuadPart) / quadpart;
    printf("notaligned_SSE_Gauss:%f\n",notaligned_SSE_Gauss);
    cout << endl;
    reset(test);

    //AVX�㷨
    QueryPerformanceCounter(&timeStart);
    avx_gauss(test);
    QueryPerformanceCounter(&timeEnd);
    double _AVX_Gauss = (timeEnd.QuadPart - timeStart.QuadPart) / quadpart;
    printf("AVX_Gauss:%f\n", _AVX_Gauss);
    cout << endl;
    reset(test);

    system("pause");
    return 0;
}
