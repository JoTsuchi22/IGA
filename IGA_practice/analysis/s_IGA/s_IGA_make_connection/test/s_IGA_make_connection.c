#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define DIMENSION 2                 // 2次元
#define MERGE_DISTANCE 1.0e-13      // コントロールポイントが同じ点と判定する距離
#define MAX_DISP_CONSTRAINT 10      // 変位指定する変位量の最大個数
#define MAX_DISP_CONSTRAINT_EDGE 10 // 変位指定する辺の最大個数
#define MAX_DISTRIBUTED_LOAD 5      // 分布荷重の最大個数

static int Total_patch;
static double E_and_nu[2];
static int disp_constraint_n[DIMENSION];
static int disp_constraint_edge_n[DIMENSION][MAX_DISP_CONSTRAINT];
static double disp_constraint_amount[DIMENSION][MAX_DISP_CONSTRAINT];
static int disp_constraint[DIMENSION][MAX_DISP_CONSTRAINT][MAX_DISP_CONSTRAINT_EDGE][3];
static int distributed_load_n;
static double distributed_load_info[MAX_DISTRIBUTED_LOAD][9];
static int counter = 0;
static int KV_to_here, CP_to_here, CP_result_to_here;
static int A_to_own, A_to_opponent, B_to_here;

void Get_inputdata_boundary(char *filename);
void Get_inputdata_patch_0(char *filename, int *temp_Order, int *temp_KV_info, int *temp_CP_info);
void Get_inputdata_patch_1(char *filename, double *temp_KV, double *temp_CP, int *temp_A, double *temp_B, int *temp_KV_info, int *temp_CP_info, int num);
void Check_B(int num_own, int num_opponent, double *temp_B, int *temp_Edge_info, int *temp_Opponent_patch_num);
void Make_connectivity(int num, int *temp_CP_info, int *temp_Edge_info, int *temp_Opponent_patch_num, int *temp_Connectivity, int *temp_A, double *temp_CP, double *temp_CP_result);
void Sort_and_Merge(int n, int *temp_A, int *temp_Boundary, int *temp_Boundary_result, int *length_before, int *length_after);
void Output_inputdata(int *temp_Order, int *temp_KV_info, int *temp_CP_info, int *temp_Connectivity, double *temp_KV, double *temp_CP_result);
// void Output_by_Gnuplot();

FILE *fp;

int main(int argc, char *argv[])
{
    int i, j, k;
    int Total_input = argc - 1;

    // ファイル読み込み
    Get_inputdata_boundary(argv[1]); // boundaryのインプットデータ処理

    // 動的メモリ確保
    int *Order = (int *)malloc(sizeof(int) * Total_patch * DIMENSION);      // int Order[パッチ番号][DIMENSION]
    int *KV_info = (int *)malloc(sizeof(int) * Total_patch * DIMENSION);    // int KV_info[パッチ番号][DIMENSION]
    int *CP_info = (int *)malloc(sizeof(int) * Total_patch * DIMENSION);    // int CP_info[パッチ番号][DIMENSION]

    if (Order == NULL || KV_info == NULL || CP_info == NULL)
    {
        printf("Memory cannot be allocated\n");
        exit(1);
    }

    // patchのインプットデータ処理(1回目)
    for (i = 1; i < Total_input; i++)
    {
        Get_inputdata_patch_0(argv[i + 1], Order, KV_info, CP_info);
    }

    int temp1 = 0; // temp1 : 全パッチ含めた総コントロールポイント数
    int temp2 = 0; // temp2 : 4辺のコントロールポイントの和を全パッチ分足した値 * 2
    int temp3 = 0; // temp3 : 全パッチ含めた総ノットベクトル数

    for (i = 0; i < Total_patch; i++)
    {
        temp1 += CP_info[i * DIMENSION] * CP_info[i * DIMENSION + 1];
        temp2 += 2 * (CP_info[i * DIMENSION] + CP_info[i * DIMENSION + 1]);
        temp3 += KV_info[i * DIMENSION] + KV_info[i * DIMENSION + 1];
    }
    printf("Total Control Point = %d\n", temp1);

    // 動的メモリ確保
    double *CP = (double *)malloc(sizeof(double) * temp1 * (DIMENSION + 1));            // double CP[パッチ番号][パッチ内CP番号][xyw -> 3]
    double *CP_result = (double *)malloc(sizeof(double) * temp1 * (DIMENSION + 1));     // double CP_result[通しのコントロールポイント番号(連番)][xyw -> 3]
    int *A = (int *)malloc(sizeof(int) * temp2);                                        // int    A[パッチ番号][辺番号(0~3)][辺内のコネクティビティ]
    double *B = (double *)malloc(sizeof(double) * Total_patch * 16 * (DIMENSION + 1));  // double B[パッチ番号][辺番号(0~7)(正負方向)][各辺の端の2頂点][座標xyw -> 3]
    int *Connectivity = (int *)malloc(sizeof(int) * temp1);                             // int    Connectivity[パッチ番号][パッチ内CP番号]
    double *KV = (double *)malloc(sizeof(double) * Total_patch * temp3);                // double KV[パッチ番号][DIMENSION][ノットベクトル番号]

    if (CP == NULL || CP_result == NULL || A == NULL || B == NULL || Connectivity == NULL || KV == NULL)
    {
        printf("Memory cannot be allocated\n");
        exit(1);
    }

    // patchのインプットデータ処理(2回目)
    counter = 0;
    KV_to_here = 0, CP_to_here = 0, B_to_here = 0;
    for (i = 1; i < Total_input; i++)
    {
        Get_inputdata_patch_1(argv[i + 1], KV, CP, A, B, KV_info, CP_info, counter);
        counter++;
    }

    printf("Done get input\n");

    // 動的メモリ確保
    int *Edge_info = (int *)malloc(sizeof(int) * Total_patch * 32);         // int Edge_info[パッチ番号][own 辺番号(正固定0~3)][opp 辺番号(0~7)]
    int *Opponent_patch_num = (int *)malloc(sizeof(int) * Total_patch * 4); // int Opponent_patch_num[パッチ番号][own 辺番号(正固定0~3]

    if (Edge_info == NULL || Opponent_patch_num == NULL)
    {
        printf("Memory cannot be allocated\n");
        exit(1);
    }

    for (i = 0; i < Total_patch * 32; i++)
    {
        Edge_info[i] = 0;
    }

    // パッチコネクティビティの作成
    printf("state: patch connectivity\n");
    counter = 0;
    CP_to_here = 0, CP_result_to_here = 0, B_to_here = 0;
    for (i = 0; i < Total_patch; i++)
    {
        for (j = 0; j < i; j++)
        {
            Check_B(i, j, B, Edge_info, Opponent_patch_num);
        }
        Make_connectivity(i, CP_info, Edge_info, Opponent_patch_num, Connectivity, A, CP, CP_result);
    }

    // 動的メモリ確保
    int temp4 = 0, temp5 = 0;
    for (i = 0; i < DIMENSION; i++)
    {
        for (j = 0; j < disp_constraint_n[i]; j++)
        {
            for (k = 0; k < disp_constraint_edge_n[i][j]; k++)
            {
                if (disp_constraint[i][j][k][1] == 0)
                {
                    temp4 += CP_info[disp_constraint[i][j][k][0] * DIMENSION + 1];
                }
                else if (disp_constraint[i][j][k][1] == 1)
                {
                    temp4 += CP_info[disp_constraint[i][j][k][0] * DIMENSION];
                }
            }
        }
        temp5 += disp_constraint_n[i];
    }

    int *length_before = (int *)malloc(sizeof(int) * temp5);    // 各変位量でのマージ前の長さ
    int *length_after = (int *)malloc(sizeof(int) * temp5);     // 各変位量でのマージ後の長さ
    int *Boundary = (int *)malloc(sizeof(int) * temp4);         // 境界条件のコネクティビティ
    int *Boundary_result = (int *)malloc(sizeof(int) * temp4);  // ソート・マージ後境界条件のコネクティビティ

    if (length_before == NULL || length_after == NULL || Boundary == NULL || Boundary_result == NULL)
    {
        printf("Memory cannot be allocated\n");
        exit(1);
    }

    // 強制変位・変位固定の境界条件を作成
    // Sort_and_Merge(temp5, A, Boundary, Boundary_result, length_before, length_after);

    printf("state: output\n");
    Output_inputdata(Order, KV_info, CP_info, Connectivity, KV, CP_result);

    // 図の出力
    // Output_by_Gnuplot();

    // メモリ解放
    free(Order), free(KV_info), free(CP_info);
    free(CP), free(CP_result), free(A), free(B), free(Connectivity), free(KV);
    free(Edge_info), free(Opponent_patch_num);
    free(length_before), free(length_after), free(Boundary), free(Boundary_result);

    return 0;
}


void Get_inputdata_boundary(char *filename)
{
    int i, j, k;
    char s[256];

    int temp_i;
    double temp_d;

    fp = fopen(filename, "r");

    // DIMENSION
    fscanf(fp, "%d", &temp_i);
    if (temp_i != 2)
    {
        printf("Error, at boundary input file\nDIMENSION must be 2 in this program\n");
        exit(1);
    }

    fgets(s, 256, fp);

    // パッチ数
    fscanf(fp, "%d", &temp_i);
    Total_patch = temp_i;
    printf("Total patch = %d\n", Total_patch);

    fgets(s, 256, fp);

    // ヤング率, ポアソン比
    for (i = 0; i < DIMENSION; i++)
    {
        fscanf(fp, "%lf", &temp_d);
        E_and_nu[i] = temp_d;
        printf("E_and_nu[%d] = %le\n", i, E_and_nu[i]);
    }

    fgets(s, 256, fp);

    // x, y方向への変位指定する個数
    for (i = 0; i < DIMENSION; i++)
    {
        fscanf(fp, "%d", &temp_i);
        disp_constraint_n[i] = temp_i;

        for (j = 0; j < disp_constraint_n[i]; j++)
        {
            fscanf(fp, "%d", &temp_i);
            disp_constraint_edge_n[i][j] = temp_i;

            fscanf(fp, "%lf", &temp_d);
            disp_constraint_amount[i][j] = temp_d;

            for (k = 0; k < disp_constraint_edge_n[i][j]; k++)
            {
                fscanf(fp, "%d", &temp_i);
                disp_constraint[i][j][k][0] = temp_i;

                fscanf(fp, "%d", &temp_i);
                disp_constraint[i][j][k][1] = temp_i;

                fscanf(fp, "%d", &temp_i);
                disp_constraint[i][j][k][2] = temp_i;
            }

            fgets(s, 256, fp);
        }
    }
    
    // 分布荷重の荷重の個数
    fscanf(fp, "%d", &temp_i);
    distributed_load_n = temp_i;

    for (i = 0; i < distributed_load_n; i++)
    {
        for (j = 0; j < 9; j++)
        {
            fscanf(fp, "%lf", &temp_d);
            distributed_load_info[i][j] = temp_d;
        }
    }

    fclose(fp);
}


void Get_inputdata_patch_0(char *filename, int *temp_Order, int *temp_KV_info, int *temp_CP_info)
{
    int i;
    char s[256];

    int temp_counter;
    int temp_i;

    printf("%s\n", filename);

    fp = fopen(filename, "r");

    // コントロールポイント数
    fscanf(fp, "%d", &temp_i);
    printf("%d\n", temp_i);

    fgets(s, 256, fp);

    // 次数
    temp_counter = counter;
    for (i = 0; i < DIMENSION; i++)
    {
        fscanf(fp, "%d", &temp_i);
        temp_Order[temp_counter] = temp_i;
        printf("temp_Order[%d] = %d\n", temp_counter, temp_Order[temp_counter]);
        temp_counter++;
    }

    fgets(s, 256, fp);

    // ノットベクトルの個数
    temp_counter = counter;
    for (i = 0; i < DIMENSION; i++)
    {
        fscanf(fp, "%d", &temp_i);
        temp_KV_info[temp_counter] = temp_i;
        printf("temp_KV_info[%d] = %d\n", temp_counter, temp_KV_info[temp_counter]);
        temp_counter++;
    }

    fgets(s, 256, fp);

    // 各方向のコントロールポイント数
    for (i = 0; i < DIMENSION; i++)
    {
        fscanf(fp, "%d", &temp_i);
        temp_CP_info[counter] = temp_i;
        printf("temp_CP_info[%d] = %d\n", counter, temp_CP_info[counter]);
        counter++;
    }

    fclose(fp);
}


void Get_inputdata_patch_1(char *filename, double *temp_KV, double *temp_CP, int *temp_A, double *temp_B, int *temp_KV_info, int *temp_CP_info, int num)
{
    int i, j;
    char s[256];

    int temp_i;
    double temp_d;

    fp = fopen(filename, "r");

    // コントロールポイント数
    fscanf(fp, "%d", &temp_i);
    printf("%d\n", temp_i);

    fgets(s, 256, fp);

    // 次数
    for (i = 0; i < DIMENSION; i++)
    {
        fscanf(fp, "%d", &temp_i);
        printf("%d\t", temp_i);
    }
    printf("\n");

    fgets(s, 256, fp);

    // ノットベクトルの個数
    for (i = 0; i < DIMENSION; i++)
    {
        fscanf(fp, "%d", &temp_i);
        printf("%d\t", temp_i);
    }
    printf("\n");

    fgets(s, 256, fp);

    // 各方向のコントロールポイント数
    for (i = 0; i < DIMENSION; i++)
    {
        fscanf(fp, "%d", &temp_i);
        printf("%d\t", temp_i);
    }
    printf("\n");

    fgets(s, 256, fp);

    // ノットベクトル
    for (i = 0; i < DIMENSION; i++)
    {
        for (j = 0; j < temp_KV_info[num * DIMENSION + i]; j++)
        {
            fscanf(fp, "%lf", &temp_d);
            temp_KV[KV_to_here + j] = temp_d;
            printf("%le", temp_KV[KV_to_here + j]);
        }
        KV_to_here += temp_KV_info[num * DIMENSION + i];
        printf("\n");
    }
    printf("\n");
    
    fgets(s, 256, fp);

    // コントロールポイント
    int temp_CP_to_here = CP_to_here;
    for (i = 0; i < temp_CP_info[num * DIMENSION] * temp_CP_info[num * DIMENSION + 1]; i++)
    {
        fscanf(fp, "%d", &temp_i);
        for (j = 0; j < DIMENSION + 1; j++)
        {
            fscanf(fp, "%lf", &temp_d);
            temp_CP[CP_to_here + j] = temp_d;
            printf("%le ", temp_CP[CP_to_here + j]);
        }
        printf("\n");
        CP_to_here += DIMENSION + 1;
    }
    printf("\n");

    printf("B_to_here = %d\n", B_to_here);

    printf("temp_CP_to_here = %d\n", temp_CP_to_here);
    printf("temp_CP_to_here = %le\n", temp_CP[temp_CP_to_here]);

    fclose(fp);
    
    // B 配列を作成
    // 辺0 点0
    temp_B[B_to_here]     = temp_CP[temp_CP_to_here];
    temp_B[B_to_here + 1] = temp_CP[temp_CP_to_here + 1];
    temp_B[B_to_here + 2] = temp_CP[temp_CP_to_here + 2];
    B_to_here += DIMENSION + 1;

    // 辺0 点1
    temp_B[B_to_here]     = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] - 1) * (DIMENSION + 1)];
    temp_B[B_to_here + 1] = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] - 1) * (DIMENSION + 1) + 1];
    temp_B[B_to_here + 2] = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] - 1) * (DIMENSION + 1) + 2];
    B_to_here += DIMENSION + 1;

    // 辺1 点0
    temp_B[B_to_here]     = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] - 1) * (DIMENSION + 1)];
    temp_B[B_to_here + 1] = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] - 1) * (DIMENSION + 1) + 1];
    temp_B[B_to_here + 2] = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] - 1) * (DIMENSION + 1) + 2];
    B_to_here += DIMENSION + 1;

    // 辺1 点1
    temp_B[B_to_here]     = temp_CP[temp_CP_to_here];
    temp_B[B_to_here + 1] = temp_CP[temp_CP_to_here + 1];
    temp_B[B_to_here + 2] = temp_CP[temp_CP_to_here + 2];
    B_to_here += DIMENSION + 1;

    // 辺2 点0
    temp_B[B_to_here]     = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] - 1) * (DIMENSION + 1)];
    temp_B[B_to_here + 1] = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] - 1) * (DIMENSION + 1) + 1];
    temp_B[B_to_here + 2] = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] - 1) * (DIMENSION + 1) + 2];
    B_to_here += DIMENSION + 1;

    // 辺2 点1
    temp_B[B_to_here]     = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] * temp_CP_info[num * DIMENSION + 1] - 1) * (DIMENSION + 1)];
    temp_B[B_to_here + 1] = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] * temp_CP_info[num * DIMENSION + 1] - 1) * (DIMENSION + 1) + 1];
    temp_B[B_to_here + 2] = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] * temp_CP_info[num * DIMENSION + 1] - 1) * (DIMENSION + 1) + 2];
    B_to_here += DIMENSION + 1;

    // 辺3 点0
    temp_B[B_to_here]     = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] * temp_CP_info[num * DIMENSION + 1] - 1) * (DIMENSION + 1)];
    temp_B[B_to_here + 1] = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] * temp_CP_info[num * DIMENSION + 1] - 1) * (DIMENSION + 1) + 1];
    temp_B[B_to_here + 2] = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] * temp_CP_info[num * DIMENSION + 1] - 1) * (DIMENSION + 1) + 2];
    B_to_here += DIMENSION + 1;

    // 辺3 点1
    temp_B[B_to_here]     = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] - 1) * (DIMENSION + 1)];
    temp_B[B_to_here + 1] = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] - 1) * (DIMENSION + 1) + 1];
    temp_B[B_to_here + 2] = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] - 1) * (DIMENSION + 1) + 2];
    B_to_here += DIMENSION + 1;

    // 辺4 点0
    temp_B[B_to_here]     = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] * (temp_CP_info[num * DIMENSION + 1] - 1)) * (DIMENSION + 1)];
    temp_B[B_to_here + 1] = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] * (temp_CP_info[num * DIMENSION + 1] - 1)) * (DIMENSION + 1) + 1];
    temp_B[B_to_here + 2] = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] * (temp_CP_info[num * DIMENSION + 1] - 1)) * (DIMENSION + 1) + 2];
    B_to_here += DIMENSION + 1;

    // 辺4 点1
    temp_B[B_to_here]     = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] * temp_CP_info[num * DIMENSION + 1] - 1) * (DIMENSION + 1)];
    temp_B[B_to_here + 1] = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] * temp_CP_info[num * DIMENSION + 1] - 1) * (DIMENSION + 1) + 1];
    temp_B[B_to_here + 2] = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] * temp_CP_info[num * DIMENSION + 1] - 1) * (DIMENSION + 1) + 2];
    B_to_here += DIMENSION + 1;

    // 辺5 点0
    temp_B[B_to_here]     = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] * temp_CP_info[num * DIMENSION + 1] - 1) * (DIMENSION + 1)];
    temp_B[B_to_here + 1] = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] * temp_CP_info[num * DIMENSION + 1] - 1) * (DIMENSION + 1) + 1];
    temp_B[B_to_here + 2] = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] * temp_CP_info[num * DIMENSION + 1] - 1) * (DIMENSION + 1) + 2];
    B_to_here += DIMENSION + 1;

    // 辺5 点1
    temp_B[B_to_here]     = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] * (temp_CP_info[num * DIMENSION + 1] - 1)) * (DIMENSION + 1)];
    temp_B[B_to_here + 1] = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] * (temp_CP_info[num * DIMENSION + 1] - 1)) * (DIMENSION + 1) + 1];
    temp_B[B_to_here + 2] = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] * (temp_CP_info[num * DIMENSION + 1] - 1)) * (DIMENSION + 1) + 2];
    B_to_here += DIMENSION + 1;

    // 辺6 点0
    temp_B[B_to_here]     = temp_CP[temp_CP_to_here];
    temp_B[B_to_here + 1] = temp_CP[temp_CP_to_here + 1];
    temp_B[B_to_here + 2] = temp_CP[temp_CP_to_here + 2];
    B_to_here += DIMENSION + 1;

    // 辺6 点1
    temp_B[B_to_here]     = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] * (temp_CP_info[num * DIMENSION + 1] - 1)) * (DIMENSION + 1)];
    temp_B[B_to_here + 1] = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] * (temp_CP_info[num * DIMENSION + 1] - 1)) * (DIMENSION + 1) + 1];
    temp_B[B_to_here + 2] = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] * (temp_CP_info[num * DIMENSION + 1] - 1)) * (DIMENSION + 1) + 2];
    B_to_here += DIMENSION + 1;

    // 辺7 点0
    temp_B[B_to_here]     = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] * (temp_CP_info[num * DIMENSION + 1] - 1)) * (DIMENSION + 1)];
    temp_B[B_to_here + 1] = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] * (temp_CP_info[num * DIMENSION + 1] - 1)) * (DIMENSION + 1) + 1];
    temp_B[B_to_here + 2] = temp_CP[temp_CP_to_here + (temp_CP_info[num * DIMENSION] * (temp_CP_info[num * DIMENSION + 1] - 1)) * (DIMENSION + 1) + 2];
    B_to_here += DIMENSION + 1;

    // 辺7 点1
    temp_B[B_to_here]     = temp_CP[temp_CP_to_here];
    temp_B[B_to_here + 1] = temp_CP[temp_CP_to_here + 1];
    temp_B[B_to_here + 2] = temp_CP[temp_CP_to_here + 2];
    B_to_here += DIMENSION + 1;


    printf("B on patch %d\n", num);
    printf("point0[x y w] point1[x y w]\n");
    int temp_counter = num * 16 * (DIMENSION + 1);

    for (i = 0; i < 8; i++)
    {
        printf("[%le %le %le]\t[%le %le %le]\n", temp_B[temp_counter], temp_B[temp_counter + 1], temp_B[temp_counter + 2], temp_B[temp_counter + 3], temp_B[temp_counter + 4], temp_B[temp_counter + 5]);
        temp_counter += 6;
    }
}


void Check_B(int num_own, int num_opponent, double *temp_B, int *temp_Edge_info, int *temp_Opponent_patch_num)
{
    int i, j;
    int ii;
    int x_diff[2], y_diff[2], w_diff[2];

    int Check_B_own_to_here = num_own * 16 * (DIMENSION + 1);
    int Check_B_opponent_to_here = num_opponent * 16 * (DIMENSION + 1);

    // printf("自分パッチの先頭x y %le %le\n", temp_B[Check_B_own_to_here], temp_B[Check_B_own_to_here + 1]);
    // printf("相手パッチの先頭x y %le %le\n", temp_B[Check_B_opponent_to_here], temp_B[Check_B_opponent_to_here + 1]);

    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 8; j++)
        {
            ii = 2 * i;
            // 点0
            x_diff[0] = temp_B[Check_B_own_to_here + ii * 2 * (DIMENSION + 1)]     - temp_B[Check_B_opponent_to_here + j * 2 * (DIMENSION + 1)];
            y_diff[0] = temp_B[Check_B_own_to_here + ii * 2 * (DIMENSION + 1) + 1] - temp_B[Check_B_opponent_to_here + j * 2 * (DIMENSION + 1) + 1];
            w_diff[0] = temp_B[Check_B_own_to_here + ii * 2 * (DIMENSION + 1) + 2] - temp_B[Check_B_opponent_to_here + j * 2 * (DIMENSION + 1) + 2];

            // 点1
            x_diff[1] = temp_B[Check_B_own_to_here + ii * 2 * (DIMENSION + 1) + (DIMENSION + 1)]     - temp_B[Check_B_opponent_to_here + j * 2 * (DIMENSION + 1) + (DIMENSION + 1)];
            y_diff[1] = temp_B[Check_B_own_to_here + ii * 2 * (DIMENSION + 1) + (DIMENSION + 1) + 1] - temp_B[Check_B_opponent_to_here + j * 2 * (DIMENSION + 1) + (DIMENSION + 1) + 1];
            w_diff[1] = temp_B[Check_B_own_to_here + ii * 2 * (DIMENSION + 1) + (DIMENSION + 1) + 2] - temp_B[Check_B_opponent_to_here + j * 2 * (DIMENSION + 1) + (DIMENSION + 1) + 2];

            // 辺が一致している場合 Edge_info を True
            if (sqrt(pow(x_diff[0], 2) + pow(y_diff[0], 2) + pow(w_diff[0], 2)) <= MERGE_DISTANCE && sqrt(pow(x_diff[1], 2) + pow(y_diff[1], 2) + pow(w_diff[1], 2)) <= MERGE_DISTANCE)
            {
                temp_Edge_info[num_own * 32 + i * 8 + j] = 1;
                temp_Opponent_patch_num[num_own * 4 + i] = num_opponent;
                printf("own_patch:%d opp_patch:%d own_edge:%d opp_edge:%d\n", num_own, num_opponent, i, j);
                return;
            }
        }
    }
}


void Make_connectivity(int num, int *temp_CP_info, int *temp_Edge_info, int *temp_Opponent_patch_num, int *temp_Connectivity, int *temp_A, double *temp_CP, double *temp_CP_result)
{
    int i, j, k;
    int p, q;
    int temp_CP_n = 0;
    int Edge[4];
    int Edge_to_here = num * 32;

    printf("make connectivity on patch %d\n", num);

    A_to_own = 0;
    for (i = 0; i < num; i++)
    {
        A_to_own += 2 * (temp_CP_info[i * DIMENSION] + temp_CP_info[i * DIMENSION + 1]);
    }

    for (i = 0; i < 4; i++)
    {
        Edge[i] = 0;
    }

    // 重なっている辺の A 配列を作成
    for (i = 0; i < 4; i++)
    {
        if (i == 0)
        {
            continue;
        }
        else if (i == 1)
        {
            A_to_own += temp_CP_info[num * DIMENSION];
        }
        else if (i == 2)
        {
            A_to_own += temp_CP_info[num * DIMENSION] + temp_CP_info[num * DIMENSION + 1];
        }
        else if (i == 3)
        {
            A_to_own += 2 * temp_CP_info[num * DIMENSION] + temp_CP_info[num * DIMENSION + 1];
        }

        for (j = 0; j < 8; j++)
        {
            if (temp_Edge_info[Edge_to_here + i * 8 + j] == 1)
            {
                A_to_opponent = 0;
                for (k = 0; k < temp_Opponent_patch_num[num * 4 + i]; k++)
                {
                    A_to_opponent += 2 * (temp_CP_info[k * DIMENSION] + temp_CP_info[k * DIMENSION + 1]);
                }

                printf("Patch num = %d\n", num);
                printf("Edge num = %d\n", j);
                printf("Opponent_patch_num = %d\n", temp_Opponent_patch_num[num * 4 + i]);

                Edge[i] = 1;

                p = j / 2;
                q = j % 2;
                if (p == 0)
                {
                    temp_CP_n = temp_CP_info[num * DIMENSION];
                }
                else if (p == 1)
                {
                    temp_CP_n = temp_CP_info[num * DIMENSION + 1];
                    A_to_opponent += temp_CP_info[temp_Opponent_patch_num[num * 4 + i] * DIMENSION];
                }
                else if (p == 2)
                {
                    temp_CP_n = temp_CP_info[num * DIMENSION];
                    A_to_opponent += temp_CP_info[temp_Opponent_patch_num[num * 4 + i] * DIMENSION] + temp_CP_info[temp_Opponent_patch_num[num * 4 + i] * DIMENSION + 1];
                }
                else if (p == 3)
                {
                    temp_CP_n = temp_CP_info[num * DIMENSION + 1];
                    A_to_opponent += 2 * temp_CP_info[temp_Opponent_patch_num[num * 4 + i] * DIMENSION] + temp_CP_info[temp_Opponent_patch_num[num * 4 + i] * DIMENSION + 1];
                }

                if (q == 0)
                {
                    for (k = 0; k < temp_CP_n; k++)
                    {
                        temp_A[A_to_own + k] = temp_A[A_to_opponent + k];
                    }
                    break;
                }
                else if (q == 1)
                {
                    for (k = 0; k < temp_CP_n; k++)
                    {
                        temp_A[A_to_own + k] = temp_A[A_to_own + (temp_CP_n - 1) - k];
                    }
                    break;
                }
            }
        }
        // Edge_to_here += 8;
    }
    
    // コネクティビティを作成
    int xi, eta;
    // int loc_num = 0;

    A_to_own = 0;
    for (i = 0; i < num; i++)
    {
        A_to_own += 2 * (temp_CP_info[i * DIMENSION] + temp_CP_info[i * DIMENSION + 1]);
    }

    for (eta = 0; eta < temp_CP_info[num * DIMENSION + 1]; eta++)
    {
        for (xi = 0; xi < temp_CP_info[num * DIMENSION]; xi++)
        {
            if (eta == 0 && Edge[0] == 1)
            {
                temp_Connectivity[CP_to_here + eta * temp_CP_info[num * DIMENSION] + xi] = temp_A[A_to_own + xi];
            }
            else if (eta == temp_CP_info[num * DIMENSION + 1] - 1 && Edge[2] == 1)
            {
                temp_Connectivity[CP_to_here + eta * temp_CP_info[num * DIMENSION] + xi] = temp_A[A_to_own + temp_CP_info[num * DIMENSION] + temp_CP_info[num * DIMENSION + 1] + xi];
            }
            else if (xi == 0 && Edge[3] == 1)
            {
                temp_Connectivity[CP_to_here + eta * temp_CP_info[num * DIMENSION] + xi] = temp_A[A_to_own + 2 * temp_CP_info[num * DIMENSION] + temp_CP_info[num * DIMENSION + 1] + eta];
            }
            else if (xi == temp_CP_info[num * DIMENSION] - 1 && Edge[1] == 1)
            {
                temp_Connectivity[CP_to_here + eta * temp_CP_info[num * DIMENSION] + xi] = temp_A[A_to_own + temp_CP_info[num * DIMENSION] + eta];
            }
            else
            {
                temp_Connectivity[CP_to_here + eta * temp_CP_info[num * DIMENSION] + xi] = counter;
                temp_CP_result[CP_result_to_here] = temp_CP[(CP_to_here + eta * temp_CP_info[num * DIMENSION] + xi) * (DIMENSION + 1)];
                temp_CP_result[CP_result_to_here + 1] = temp_CP[(CP_to_here + eta * temp_CP_info[num * DIMENSION] + xi) * (DIMENSION + 1) + 1];
                temp_CP_result[CP_result_to_here + 2] = temp_CP[(CP_to_here + eta * temp_CP_info[num * DIMENSION] + xi) * (DIMENSION + 1) + 2];
                counter++;
                CP_result_to_here += (DIMENSION + 1);
            }

            // A 配列の作ってない分を作成
            if (eta == 0 && Edge[0] == 0)
            {
                temp_A[A_to_own + xi] = temp_Connectivity[CP_to_here + eta * temp_CP_info[num * DIMENSION] + xi];
            }
            if (eta == temp_CP_info[num * DIMENSION + 1] - 1 && Edge[2] == 0)
            {
                temp_A[A_to_own + temp_CP_info[num * DIMENSION] + temp_CP_info[num * DIMENSION + 1] + xi] = temp_Connectivity[CP_to_here + eta * temp_CP_info[num * DIMENSION] + xi];
            }
            if (xi == 0 && Edge[3] == 0)
            {
                temp_A[A_to_own + 2 * temp_CP_info[num * DIMENSION] + temp_CP_info[num * DIMENSION + 1] + eta] = temp_Connectivity[CP_to_here + eta * temp_CP_info[num * DIMENSION] + xi];
            }
            if (xi == temp_CP_info[num * DIMENSION] - 1 && Edge[1] == 0)
            {
                temp_A[A_to_own + temp_CP_info[num * DIMENSION] + eta] = temp_Connectivity[CP_to_here + eta * temp_CP_info[num * DIMENSION] + xi];
            }
        }
    }
    CP_to_here += temp_CP_info[num * DIMENSION] * temp_CP_info[num * DIMENSION + 1];
}


// void Sort_and_Merge(int n, int *temp_A, int *temp_Boundary, int *temp_Boundary_result, int *length_before, int *length_after)
// {
//     int i, j, k, l, m, n;
//     int temp = 0;

//     for (i = 0; i < DIMENSION; i++)
//     {
//         for (j = 0; j < disp_constraint_n[i]; j++)
//         {
//             for (k = 0; k < disp_constraint_edge_n[i][j]; k++)
//             {
//                 if (disp_constraint[i][j][k][1] == 0)
//                 {
//                     length_before[temp] += CP_info[disp_constraint[i][j][k][0] * DIMENSION + 1];
//                 }
//                 else if (disp_constraint[i][j][k][1] == 1)
//                 {
//                     length_before[temp] += CP_info[disp_constraint[i][j][k][0] * DIMENSION];
//                 }
//             }
//             temp++;
//         }
//     }

//     temp = 0;

//     for (i = 0; i < DIMENSION; i++)
//     {
//         for (j = 0; j < disp_constraint_n[i]; j++)
//         {
//             for (k = 0; k < disp_constraint_edge_n[i][j]; k++)
//             {
//                 int A_to_here = 0;
//                 for (l = 0; l < disp_constraint[i][j][k][0]; l++)
//                 {
//                     A_to_here += 2 * (CP_info[l * DIMENSION] + CP_info[l * DIMENSION + 1]);
//                 }

//                 if (disp_constraint[i][j][k][1] == 1 || disp_constraint[i][j][k][2] == 0)
//                 {
//                     continue;
//                 }
//                 else if (disp_constraint[i][j][k][1] == 0 || disp_constraint[i][j][k][2] == 1)
//                 {
//                     A_to_here += CP_info[disp_constraint[i][j][k][0] * DIMENSION];
//                 }
//                 else if (disp_constraint[i][j][k][1] == 1 || disp_constraint[i][j][k][2] == 1)
//                 {
//                     A_to_here += CP_info[disp_constraint[i][j][k][0] * DIMENSION] + CP_info[disp_constraint[i][j][k][0] * DIMENSION + 1];
//                 }
//                 else if (disp_constraint[i][j][k][1] == 0 || disp_constraint[i][j][k][2] == 0)
//                 {
//                     A_to_here += 2 * CP_info[disp_constraint[i][j][k][0] * DIMENSION] + CP_info[disp_constraint[i][j][k][0] * DIMENSION + 1];
//                 }
                
//                 if (disp_constraint[i][j][k][1] == 0)
//                 {
//                     for (l = 0; l < CP_info[disp_constraint[i][j][k][0] * DIMENSION] + 1; l++)
//                     {
//                         temp_Boundary[temp] = temp_A[A_to_here + l];
//                         temp++;
//                     }
//                 }
//                 else if (disp_constraint[i][j][k][1] == 1)
//                 {
//                     for (l = 0; l < CP_info[disp_constraint[i][j][k][0] * DIMENSION]; l++)
//                     {
//                         temp_Boundary[temp] = temp_A[A_to_here + l];
//                         temp++;
//                     }
//                 }
//             }
//         }
//     }


// }


void Output_inputdata(int *temp_Order, int *temp_KV_info, int *temp_CP_info, int *temp_Connectivity, double *temp_KV, double *temp_CP_result)
{
    int i, j, k;
    char str[256] = "input.txt";
    
    fp = fopen(str, "w");

    // ヤング率
    fprintf(fp, "%d\t", (int)E_and_nu[0]);

    // ポアソン比
    fprintf(fp, "%le\n\n", E_and_nu[1]);

    // パッチ数
    fprintf(fp, "%d\n\n", Total_patch);

    // コントロールポイント数
    fprintf(fp, "%d\n\n", (CP_result_to_here + 1) / 3);

    // 各パッチ内での各方向の次数
    for (i = 0; i < Total_patch; i++)
    {
        for (j = 0; j < DIMENSION; j++)
        {
            if (j == 0)
            {
                fprintf(fp, "%d", temp_Order[i * DIMENSION + j]);
            }
            else
            {
                fprintf(fp, "\t%d", temp_Order[i * DIMENSION + j]);
            }
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // 各パッチ内での各方向のノットベクトルの数
    for (i = 0; i < Total_patch; i++)
    {
        for (j = 0; j < DIMENSION; j++)
        {
            if (j == 0)
            {
                fprintf(fp, "%d", temp_KV_info[i * DIMENSION + j]);
            }
            else
            {
                fprintf(fp, "\t%d", temp_KV_info[i * DIMENSION + j]);
            }
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // 各パッチ内での各方向のコントロールポイントの数
    for (i = 0; i < Total_patch; i++)
    {
        for (j = 0; j < DIMENSION; j++)
        {
            if (j == 0)
            {
                fprintf(fp, "%d", temp_CP_info[i * DIMENSION + j]);
            }
            else
            {
                fprintf(fp, "\t%d", temp_CP_info[i * DIMENSION + j]);
            }
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // パッチコネクティビティ
    CP_to_here = 0;
    for (i = 0; i < Total_patch; i++)
    {
        for (j = 0; j < temp_CP_info[i * DIMENSION] * temp_CP_info[i * DIMENSION + 1]; j++)
        {
            if (j == 0)
            {
                fprintf(fp, "%d", temp_Connectivity[CP_to_here + j]);
            }
            else
            {
                fprintf(fp, "\t%d", temp_Connectivity[CP_to_here + j]);
            }
        }
        CP_to_here += temp_CP_info[i * DIMENSION] * temp_CP_info[i * DIMENSION + 1];
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // 変位拘束するコントロールポイントの数
    // fprintf(fp, "%d", XX);
    //

    // 荷重条件を与えるコントロールポイントの数
    fprintf(fp, "\t%d", 0);

    // 分布荷重の数
    fprintf(fp, "\t%d\n\n", distributed_load_n);

    // 各パッチでの各方向のノットベクトル
    KV_to_here = 0;
    for (i = 0; i < Total_patch; i++)
    {
        for (j = 0; j < DIMENSION; j++)
        {
            for (k = 0; k < temp_KV_info[i * DIMENSION + j]; k++)
            {
                if (k == 0)
                {
                    fprintf(fp, "%.16e", temp_KV[KV_to_here + k]);
                }
                else
                {
                    fprintf(fp, "\t%.16e", temp_KV[KV_to_here + k]);
                }
            }
            KV_to_here += temp_KV_info[i * DIMENSION + j];
            fprintf(fp, "\n");
        }
    }
    fprintf(fp, "\n");

    // コントロールポイント
    for (i = 0; i < (CP_result_to_here + 1) / 3; i++)
    {
        fprintf(fp, "%d", i);
        fprintf(fp, "\t%.16e", temp_CP_result[i * 3]);
        fprintf(fp, "\t%.16e", temp_CP_result[i * 3 + 1]);
        fprintf(fp, "\t%.16e\n", temp_CP_result[i * 3 + 2]);
    }
    fprintf(fp, "\n");

    // 拘束するコントロールポイント
    // for (i = 0; i < DIMENSION; i++)
    // {
    //     for (j = 0; j < disp_constraint_n[i])
    //     {
    //         for (k = 0; k < length_before[]; k++)
    //         {
    //             fprintf(fp, "%d", temp_Boundary[boundary_to_here + ]);
    //             fprintf(fp, "\t%d", i);
    //             fprintf(fp, "\t%.16e\n", disp_constraint_amount[i][j]);
    //         }
    //         boundary_to_here += length[i][j];
    //     }
    // }

    // 分布荷重
    for (i = 0; i < distributed_load_n; i++)
    {
        if (i != 0)
        {
            fprintf(fp, "\n");
        }
        fprintf(fp, "%d", (int)distributed_load_info[i][0]);
        fprintf(fp, "\t%d", (int)distributed_load_info[i][1]);
        fprintf(fp, "\t%d", (int)distributed_load_info[i][2]);
        fprintf(fp, "\t%le", distributed_load_info[i][3]);
        fprintf(fp, "\t%le", distributed_load_info[i][4]);
        fprintf(fp, "\t%le", distributed_load_info[i][5]);
        fprintf(fp, "\t%le", distributed_load_info[i][6]);
        fprintf(fp, "\t%le", distributed_load_info[i][7]);
        fprintf(fp, "\t%le", distributed_load_info[i][8]);
    }

    fclose(fp);
}


// void Output_by_Gnuplot()
// {

// }