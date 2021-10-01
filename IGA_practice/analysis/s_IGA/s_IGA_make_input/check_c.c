#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAX_N_INPUTFILE 5
#define MAX_DIMENSION 3
#define MAX_ORDER 10
#define MAX_N_Controlpoint_in_Patch 10000
#define MAX_N_Controlpoint_each_parameter 100
#define MAX_N_KNOT 1000

static int Dimension[MAX_N_INPUTFILE];
static int Total_Control_Point[MAX_N_INPUTFILE];
static int Order[MAX_N_INPUTFILE][MAX_DIMENSION];
static int knot_n[MAX_N_INPUTFILE][MAX_DIMENSION];
static int Control_point_n[MAX_N_INPUTFILE][MAX_DIMENSION];
static double knot[MAX_N_INPUTFILE][MAX_DIMENSION][MAX_N_KNOT];
static double x[MAX_N_INPUTFILE][MAX_N_Controlpoint_in_Patch];
static double y[MAX_N_INPUTFILE][MAX_N_Controlpoint_in_Patch];
static double w[MAX_N_INPUTFILE][MAX_N_Controlpoint_in_Patch];
static int OE_n[MAX_N_INPUTFILE][MAX_DIMENSION];
static int KI_uniform_n[MAX_N_INPUTFILE][MAX_DIMENSION];
static int KI_non_uniform_n[MAX_N_INPUTFILE][MAX_DIMENSION];
static double insert_knot[MAX_N_INPUTFILE][MAX_DIMENSION][MAX_N_KNOT];
static double insert_knot_in_KI[MAX_N_INPUTFILE][MAX_DIMENSION][MAX_N_KNOT];
static double removal_knot[MAX_N_INPUTFILE][MAX_DIMENSION][MAX_N_KNOT];
static int vec_length1[MAX_N_INPUTFILE][MAX_DIMENSION];
static int vec_length2[MAX_N_INPUTFILE][MAX_DIMENSION];
static double temp_knot1[MAX_N_KNOT];
static double temp_knot2[MAX_N_KNOT];
static double x_array[MAX_N_Controlpoint_each_parameter][MAX_N_Controlpoint_each_parameter];
static double y_array[MAX_N_Controlpoint_each_parameter][MAX_N_Controlpoint_each_parameter];
static double w_array[MAX_N_Controlpoint_each_parameter][MAX_N_Controlpoint_each_parameter];
static double temp_x_array[MAX_N_Controlpoint_each_parameter];
static double temp_y_array[MAX_N_Controlpoint_each_parameter];
static double temp_w_array[MAX_N_Controlpoint_each_parameter];
static double temp_x1[MAX_N_Controlpoint_in_Patch];
static double temp_y1[MAX_N_Controlpoint_in_Patch];
static double temp_w1[MAX_N_Controlpoint_in_Patch];
static double temp_x2[MAX_N_Controlpoint_in_Patch];
static double temp_y2[MAX_N_Controlpoint_in_Patch];
static double temp_w2[MAX_N_Controlpoint_in_Patch];
static int temp_Bezier_array[MAX_N_Controlpoint_each_parameter];
static int counter;
static double Bezier_x[MAX_N_Controlpoint_each_parameter][MAX_ORDER];
static double Bezier_y[MAX_N_Controlpoint_each_parameter][MAX_ORDER];
static double Bezier_w[MAX_N_Controlpoint_each_parameter][MAX_ORDER];
static int number_of_Bezier_line;

void Calc_insert_knot_in_OE(int tm, int elevation_parameter_axis)
{
    int m, n;

    m = knot_n[tm][elevation_parameter_axis];
    n = Order[tm][elevation_parameter_axis];

    double vec_temp1[MAX_N_KNOT];
    double vec_temp2[MAX_N_KNOT];
    double vec_unique[MAX_N_KNOT];

    int i, j, k;
    int temp1 = 0;
    for (i = 1; i < m - 1; i++)
    {
        if (knot[tm][elevation_parameter_axis][i - 1] != knot[tm][elevation_parameter_axis][i] && knot[tm][elevation_parameter_axis][i] != knot[tm][elevation_parameter_axis][i + 1])
        {
            vec_unique[temp1] = knot[tm][elevation_parameter_axis][i];
            temp1++;
        }
    }

    int temp3 = 0;
    for (i = 0; i < temp1; i++)
    {
        int temp2 = 0;
        for (j = 0; j < m; j++)
        {
            if (vec_unique[i] == knot[tm][elevation_parameter_axis][j])
            {
                temp2++;
            }
        }
        if (temp2 < n)
        {
            vec_temp1[temp3] = vec_unique[i];
            temp3++;
        }
    }

    int temp4 = 0;
    if (n - 2 > 0)
    {
        for (j = 0; j < temp3; j++)
        {
            for (k = 0; k < n -1; k++)
            {
                vec_temp2[temp4] = vec_temp1[j];
                temp4++;
            }
        }
        printf("temp4 = %d", temp4);
    }

    if (n - 2 > 0)
    {
        for (i = 0; i < temp4; i++)
        {
            insert_knot_in_KI[tm][elevation_parameter_axis][i] = vec_temp2[i];
        }
        vec_length1[tm][elevation_parameter_axis] = temp4;
    }
    else
    {
        for (i = 0; i < temp3; i++)
        {
            insert_knot_in_KI[tm][elevation_parameter_axis][i] = vec_temp1[i];
        }
        vec_length1[tm][elevation_parameter_axis] = temp3;
    }

    for (i = 0; i < vec_length1[tm][elevation_parameter_axis]; i++)
    {
        printf("%lf\t", insert_knot_in_KI[tm][elevation_parameter_axis][i]);
    }
    printf("\n");

    for (i = 0; i < temp3; i++)
    {
        removal_knot[tm][elevation_parameter_axis][i] = vec_temp1[i];
    }
    vec_length2[tm][elevation_parameter_axis] = temp3;
}

void KI_calc_knot_1D(int tm, int insert_parameter_axis)
{
    int i;
    for (i = 0; i < knot_n[tm][insert_parameter_axis]; i++)
    {
        temp_knot1[i] = knot[tm][insert_parameter_axis][i];
    }

    int temp1 = 0, temp2 = 0;
    for (i = 0; i < knot_n[tm][insert_parameter_axis] - 1; i++)
    {
        temp_knot2[temp1] = temp_knot1[i];
        temp1++;
        if (temp_knot1[i] < insert_knot_in_KI[tm][insert_parameter_axis][temp2] && insert_knot_in_KI[tm][insert_parameter_axis][temp2] <= temp_knot1[i + 1])
        {
            temp_knot2[temp1] = insert_knot_in_KI[tm][insert_parameter_axis][temp2];
            temp1++;
            temp2++;
        }
        if (i == knot_n[tm][insert_parameter_axis] - 2)
        {
            temp_knot2[temp1] = temp_knot1[i+1];
            temp1++;
        }
    }
}


int main()
{
    // knot[0][0][0] = 0.0;
    // knot[0][0][1] = 0.0;
    // knot[0][0][2] = 0.0;
    // knot[0][0][3] = 0.2;
    // knot[0][0][4] = 0.4;
    // knot[0][0][5] = 0.6;
    // knot[0][0][6] = 0.8;
    // knot[0][0][7] = 1.0;
    // knot[0][0][8] = 1.0;
    // knot[0][0][9] = 1.0;

    // knot[0][1][0] = 0.0;
    // knot[0][1][1] = 0.0;
    // knot[0][1][2] = 0.0;
    // knot[0][1][3] = 0.2;
    // knot[0][1][4] = 0.4;
    // knot[0][1][5] = 0.6;
    // knot[0][1][6] = 0.8;
    // knot[0][1][7] = 1.0;
    // knot[0][1][8] = 1.0;
    // knot[0][1][9] = 1.0;

    // Order[0][0] = 2;
    // Order[0][1] = 2;

    // knot_n[0][0] = 10;
    // knot_n[0][1] = 10;


    // knot[0][0][0] = 0.0;
    // knot[0][0][1] = 0.0;
    // knot[0][0][2] = 0.5;
    // knot[0][0][3] = 1.0;
    // knot[0][0][4] = 1.0;

    // knot[0][1][0] = 0.0;
    // knot[0][1][1] = 0.0;
    // knot[0][1][2] = 0.5;
    // knot[0][1][3] = 1.0;
    // knot[0][1][4] = 1.0;

    // Order[0][0] = 1;
    // Order[0][1] = 1;

    // knot_n[0][0] = 5;
    // knot_n[0][1] = 5;


    knot[0][0][0] = 0.0;
    knot[0][0][1] = 0.0;
    knot[0][0][2] = 0.0;
    knot[0][0][3] = 0.0;
    knot[0][0][4] = 0.25;
    knot[0][0][5] = 0.5;
    knot[0][0][6] = 0.75;
    knot[0][0][7] = 1.0;
    knot[0][0][8] = 1.0;
    knot[0][0][9] = 1.0;
    knot[0][0][10] = 1.0;

    knot[0][1][0] = 0.0;
    knot[0][1][1] = 0.0;
    knot[0][1][2] = 0.0;
    knot[0][1][3] = 0.0;
    knot[0][1][4] = 0.25;
    knot[0][1][5] = 0.5;
    knot[0][1][6] = 0.75;
    knot[0][1][7] = 1.0;
    knot[0][1][8] = 1.0;
    knot[0][1][9] = 1.0;
    knot[0][1][10] = 1.0;

    Order[0][0] = 3;
    Order[0][1] = 3;

    knot_n[0][0] = 11;
    knot_n[0][1] = 11;

    Calc_insert_knot_in_OE(0, 0);
    Calc_insert_knot_in_OE(0, 1);

    // int i, j;
    // for (i = 0; i < 2; i++)
    // {
    //     for (j = 0; j < 100; j++)
    //     {
    //         printf("%d\t%d\t%le\n", i,j,insert_knot[0][i][j]);
    //     }
    // }

    return 0;
}