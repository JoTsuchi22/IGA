/*****************************************************************************
 * s_IGA_make_input_crack.c :
 *****************************************************************************
 * s_IGAにおける，き裂の特異パッチのインプットデータを作成するためのプログラムです
 * このプログラムのインプットデータについてはREADMEを読んでください :)
 * 
 * 2022.4.14
 *****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define DIMENSION 2                 // 2次元
// #define MERGE_DISTANCE 1.0e-13      // コントロールポイントが同じ点と判定する距離
// #define MAX_DISP_CONSTRAINT 10      // 変位指定する変位量の最大個数
// #define MAX_DISP_CONSTRAINT_EDGE 10 // 変位指定する辺の最大個数
// #define MAX_DISTRIBUTED_LOAD 5      // 分布荷重の最大個数

FILE *fp;

void Get_input_glo(char *filename);
void Get_input_loc(char *filename);



int main(int argc, char **argv)
{
    // ファイル読み込み
    Get_input_glo(argv[1]);
    Get_input_loc(argv[2]);

    int i, j;

    for (i = 0; i < argc - 1; i++)
    {
        printf("df\n");
    }



    return 0;
}


void Get_input_glo(char *filename)
{
    fp = fopen(filename, "r");
    fclose(fp);
}


void Get_input_loc(char *filename)
{
    fp = fopen(filename, "r");
    fclose(fp);
}