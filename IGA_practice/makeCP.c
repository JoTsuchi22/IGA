#include <stdio.h>
int main(void){
    int k;
    float a = 1.0;
    int i,j;
    float x[12] = {0.00,5.00,15.00,25.00,35.00,45.00,55.00,65.00,75.00,85.00,95.00,100.00}, y[6] = {0.00,1.25,3.75,6.25,8.75,10.00};
    for (i = 0;i < 6;i++){
        for (j = 0;j < 12;j++){
            k = j + 12 * i;
            printf("%d  %f  %f  %f\n",
                k,
                x[j],
                y[i],
                a);
        }
    }
}