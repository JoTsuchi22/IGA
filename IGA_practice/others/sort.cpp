#include <stdio.h>
#include <iostream>
#include <string>

#define n_MAX_NUM 10000
#define m_MAX_NUM 4

double numberelements[n_MAX_NUM][m_MAX_NUM];

using namespace std;
int Input(char *filename)
{
    FILE *fp;
    int i;

    cout << "Reading data from " << filename << "." << endl;

    fp = fopen(filename, "r");

    if (fp == NULL) {
        cout << "Error: Input file " << filename << " is not found." << endl;
        exit(1);
    }
    
    for (i = 0; i < n_MAX_NUM; i++) {
        fscanf(fp, "%lf %lf %lf %lf",
            &numberelements[i][0],
            &numberelements[i][1],
            &numberelements[i][2],
            &numberelements[i][3]);
    }

    fclose(fp);

    return 0;
}

int Output1(char *filename)
{
    FILE *fp;
    int i;

    cout << "Wrinting data to " << filename << "." << endl;

    fp = fopen(filename, "w");    
    
    for (i = 0; i < n_MAX_NUM; i++) {
        if (numberelements[i][1] == 0){
            if (i != 0 && numberelements[i][0] == 0 && numberelements[i][1] == 0 && numberelements[i][2] == 0 && numberelements[i][3] == 0){
                continue;
            }
            else{
                fprintf(fp, "%f %f  %f  %f\n",
                    numberelements[i][0],
                    numberelements[i][1],
                    numberelements[i][2],
                    numberelements[i][3]);
            }
        }
    }
    
    fclose(fp);

    return 0;
}

int Output2(char *filename)
{
    FILE *fp;
    int i;

    cout << "Wrinting data to " << filename << "." << endl;

    fp = fopen(filename, "w");    
    
    for (i = 0; i < n_MAX_NUM; i++) {
        if (numberelements[i][1] == 5){
            fprintf(fp, "%f %f  %f  %f\n",
                numberelements[i][0],
                numberelements[i][1],
                numberelements[i][2],
                numberelements[i][3]);
        }
    }
    
    fclose(fp);

    return 0;
}

int Output3(char *filename)
{
    FILE *fp;
    int i;

    cout << "Wrinting data to " << filename << "." << endl;

    fp = fopen(filename, "w");    
    
    for (i = 0; i < n_MAX_NUM; i++) {
        if (numberelements[i][1] == 10){
            fprintf(fp, "%f %f  %f  %f\n",
                numberelements[i][0],
                numberelements[i][1],
                numberelements[i][2],
                numberelements[i][3]);
        }
    }
    
    fclose(fp);

    return 0;
}

int main(int argc, char *argv[])
{
    Input(argv[1]);
    Output1(argv[2]);
    Output2(argv[3]);
    Output3(argv[4]);
    return 0;
}