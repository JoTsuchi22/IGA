#include <stdio.h>
// #include <string.h>
// #include <math.h>
#include <stdlib.h>

#define MERGE_DISTANCE 1.0e-13 // コントロールポイントが同じ点と判定する距離

void Get_inputdata(int tm, char *filename);
void 
//図の出力
void Output_by_Gnuplot();
void Output_inputdata();

FILE *fp;

int main(int argc, char *argv[])
{
    Get_inputdata();

    // 動的メモリ確保

    Output_by_Gnuplot();

    Output_inputdata();

    // メモリ解放
    free(), free()

}

void Get_inputdata(int tm, char *filename)
{

}


void 
{

}


void Output_by_Gnuplot()
{

}


void Output_inputdata()
{
    
}