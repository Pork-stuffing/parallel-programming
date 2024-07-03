//# include <sys/time.h>
#include <windows.h>
#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;


//
//unsigned int Act[8399][264] = { 0 };
//unsigned int Pas[8399][264] = { 0 };
//
//const int Num = 263;
//const int pasNum = 4535;
//const int lieNum = 8399;




//
//unsigned int Act[23045][722] = { 0 };
//unsigned int Pas[23045][722] = { 0 };
//
//const int Num = 721;
//const int pasNum = 14325;
//const int lieNum = 23045;



unsigned int Act[37960][1188] = { 0 };
unsigned int Pas[37960][1188] = { 0 };

const int Num = 1187;
const int pasNum = 14921;
const int lieNum = 37960;



/*
unsigned int Act[43577][1363] = { 0 };
unsigned int Pas[54274][1363] = { 0 };

const int Num = 1362;
const int pasNum = 54274;
const int lieNum = 43577;
*/


//��Ԫ�ӳ�ʼ��
void init_A()
{
    //ÿ����Ԫ�ӵ�һ��Ϊ1λ���ڵ�λ�ã����������ڶ�ά������к�
    //���磺��Ԫ�ӣ�561��...����Act[561][]���
    unsigned int a;
    ifstream infile("act2.txt");
    char fin[10000] = { 0 };
    int index;
    //���ļ�����ȡ��
    while (infile.getline(fin, sizeof(fin)))
    {
        std::stringstream line(fin);
        int biaoji = 0;

        //��������ȡ����������
        while (line >> a)
        {
            if (biaoji == 0)
            {
                //ȡÿ�е�һ������Ϊ�б�
                index = a;
                biaoji = 1;
            }
            int k = a % 32;
            int j = a / 32;

            int temp = 1 << k;
            Act[index][Num - 1 - j] += temp;
            Act[index][Num] = 1;//���ø�λ�ü�¼��Ԫ�Ӹ����Ƿ�Ϊ�գ�Ϊ������0������Ϊ1
        }
    }
}

//����Ԫ�г�ʼ��
void init_P()
{
    //ֱ�Ӱ��մ����ļ���˳��棬�ڴ����ļ��ǵڼ��У���������ǵڼ���
    unsigned int a;
    ifstream infile("pas2.txt");
    char fin[10000] = { 0 };
    int index = 0;
    //���ļ�����ȡ��
    while (infile.getline(fin, sizeof(fin)))
    {
        std::stringstream line(fin);
        int biaoji = 0;

        //��������ȡ����������
        while (line >> a)
        {
            if (biaoji == 0)
            {
                //��Pas[ ][263]��ű���Ԫ��ÿ�е�һ�����֣�����֮�����Ԫ����
                Pas[index][Num] = a;
                biaoji = 1;
            }

            int k = a % 32;
            int j = a / 32;

            int temp = 1 << k;
            Pas[index][Num - 1 - j] += temp;
        }
        index++;
    }
}


void f_ordinary()
{
    int i;
    for (i = lieNum - 1; i - 8 >= -1; i -= 8)
    {
        //ÿ�ִ���8����Ԫ�ӣ���Χ�������� i-7 �� i

        for (int j = 0; j < pasNum; j++)
        {
            //��4535������Ԫ����û�������ڴ˷�Χ�ڵ�
            while (Pas[j][Num] <= i && Pas[j][Num] >= i - 7)
            {
                int index = Pas[j][Num];
                if (Act[index][Num] == 1)//��Ԫ�Ӳ�Ϊ��
                {
                    //Pas[j][]��Act[��Pas[j][18]��][]�����
                    for (int k = 0; k < Num; k++)
                    {
                        Pas[j][k] = Pas[j][k] ^ Act[index][k];
                    }

                    //����Pas[j][18]�������ֵ
                    //�������֮������������������浽Pas[j][18]�������ڷ�Χ������whileѭ��
                    //�����֮��Pas[j][ ]������
                    int num = 0, S_num = 0;
                    for (num = 0; num < Num; num++)
                    {
                        if (Pas[j][num] != 0)
                        {
                            unsigned int temp = Pas[j][num];
                            while (temp != 0)
                            {
                                temp = temp >> 1;
                                S_num++;
                            }
                            S_num += num * 32;
                            break;
                        }
                    }
                    Pas[j][Num] = S_num - 1;

                }
                else//��Ԫ��Ϊ��
                {
                    //Pas[j][]��������Ԫ��
                    for (int k = 0; k < Num; k++)
                        Act[index][k] = Pas[j][k];

                    Act[index][Num] = 1;//������Ԫ�ӷǿ�
                    break;
                }

            }
        }
    }


    for (int i = lieNum % 8 - 1; i >= 0; i--)
    {
        //ÿ�ִ���1����Ԫ�ӣ���Χ���������i

        for (int j = 0; j < pasNum; j++)
        {
            //��53������Ԫ����û���������i��
            while (Pas[j][Num] == i)
            {
                if (Act[i][Num] == 1)//��Ԫ�Ӳ�Ϊ��
                {
                    //Pas[j][]��Act[i][]�����
                    for (int k = 0; k < Num; k++)
                    {
                        Pas[j][k] = Pas[j][k] ^ Act[i][k];
                    }

                    //����Pas[j][18]�������ֵ
                    //�������֮������������������浽Pas[j][18]�������ڷ�Χ������whileѭ��
                    //�����֮��Pas[j][ ]������
                    int num = 0, S_num = 0;
                    for (num = 0; num < Num; num++)
                    {
                        if (Pas[j][num] != 0)
                        {
                            unsigned int temp = Pas[j][num];
                            while (temp != 0)
                            {
                                temp = temp >> 1;
                                S_num++;
                            }
                            S_num += num * 32;
                            break;
                        }
                    }
                    Pas[j][Num] = S_num - 1;

                }
                else//��Ԫ��Ϊ��
                {
                    //Pas[j][]��������Ԫ��
                    for (int k = 0; k < Num; k++)
                        Act[i][k] = Pas[j][k];

                    Act[i][Num] = 1;//������Ԫ�ӷǿ�
                    break;
                }
            }
        }
    }
}

int main()
{
    
    //struct timeval head, tail;
    double seconds;
    long long head, tail, freq, noww;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);


    init_A();
    init_P();
    //gettimeofday(&head, NULL);//��ʼ��ʱ
    QueryPerformanceCounter((LARGE_INTEGER*)&head);//��ʼ��ʱ
    f_ordinary();
    //gettimeofday(&tail, NULL);//������ʱ
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);//������ʱ
    seconds = (tail - head) * 1000.0 / freq;//��λ ms
    cout << seconds ;
    cout << "ms";

    //double seconds2 = ((tail.tv_sec - head.tv_sec) * 1000000 + (tail.tv_usec - head.tv_usec)) / 1000.0;//��λ ms
    //cout << "f_ordinary: " << seconds2 << " ms" << endl;


}





