/***************************************************************************************

rhinocerosのigsファイルからigaのinputデータを作るコード(マルチパッチ対応)

[argument]
argv[1]・・・igsファイル（Rhinocerosのoutputデータ）
argv[2]・・・txtファイル（今回のoutputデータのファイル名＝IGAのInputデータ）

解析をするときは
このコードで作られるtxtファイルの
    パッチ番号
    変位・集中荷重・分布荷重拘束の数
    変位・集中荷重・分布荷重拘束
を自分で入力してから行う

***************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#define ERROR					-999
#define MAX_PATCH               6

static int number1_a;
static int P_next;
static int time_next[100];
static int no_128[MAX_PATCH];
static int patch_no;
static int No_con_xi_1[MAX_PATCH],No_con_eta_1[MAX_PATCH],No_con_xi[MAX_PATCH],No_con_eta[MAX_PATCH];
static int No_knot_xi[MAX_PATCH],No_knot_eta[MAX_PATCH];
static int Order_xi[MAX_PATCH],Order_eta[MAX_PATCH];
static int all_number[MAX_PATCH];
static double D_number[MAX_PATCH][1000];
static double knot_xi[MAX_PATCH][1000],knot_eta[MAX_PATCH][1000];
static double weight[MAX_PATCH][1000];
static double con_coord[MAX_PATCH][1000][1000];

void read_igs(char input_file_name[])
{
    FILE *input_file;
    input_file=fopen(input_file_name,"r");
    assert(input_file != NULL);

    char temp_char[300],*p;
    int j,jjjj;
    int k;
    int i,ii,iii,iiii;
    int l,ll;
    int n,nn;
    int m;
    int noP;
    int number5_a,number5_b,number5_c,number5_d,number5_e;
    double number;
    int power;
    char D;

    char time[]="10,0,13H";
    char D0D[]="0D";
    char char_P[9];

    while(fgets(temp_char,300,input_file) != NULL)
    {
        p = strstr(temp_char,time);
        if(p != NULL)
        {
            nn=0;
            fgets(temp_char,200,input_file);
            fscanf(input_file,"%d",&time_next[nn]);
            printf("time_next[%d]=%d\n",nn,time_next[nn]);
            fgets(temp_char,200,input_file);
            //fgets(temp_char,200,input_file);
            nn++;
            m=0;
            for(n=nn;n<50;n++){
                p=strstr(temp_char,D0D);
                if(p!=NULL){
                    fgets(temp_char,200,input_file);
                    fscanf(input_file,"%d",&time_next[n]);
                    printf("time_next[%d]=%d\n",n,time_next[n]);
                    if(time_next[n]==128){
                        no_128[m]=n;
                        m++;
                    }
                    else{
                        m=m;
                    }
                    printf("m=%d\n",m);
                    fgets(temp_char,200,input_file);
                    //fgets(temp_char,200,input_file);
                }
                else{
                    break;
                }
            } 
            printf("time_next had read![%d]\n",n-2);
            //printf("patch_no=%d\n",m);
            patch_no=m;
            printf("patch_no=%d\n",patch_no);
            //break;
            for(l=0;l<patch_no;l++){
                //printf("no_128[%d]=%d\n",l,no_128[l]);
                noP=no_128[l]*2+1;
                //printf("noP=%d\n",noP);
                //printf("char_P=%7.dP\n",noP);
                sprintf(char_P,"00000%dP",noP);
                printf("char_P=%s\n",char_P);
                ll=0;
                while(fgets(temp_char,300,input_file) != NULL){
                    p=strstr(temp_char,char_P);
                    if(p != NULL){
                        fgets(temp_char,200,input_file);
                        //printf("%s",temp_char);
                        fscanf(input_file,"%d\t",&P_next);
                        //printf("P_next=%d\n",number1_a);
                        fseek(input_file,-170,SEEK_CUR);
                        fgets(temp_char,200,input_file);
                        //printf("%s",temp_char);
                        fscanf(input_file,"%d,",&number1_a);
                        printf("number1_a=%d\n",number1_a);
                        //printf("ll=%d\n",ll);
                        if(number1_a != 128){
                            if(ll+1==patch_no){
                                printf("Finish!!\n");
                            } 
                            else{
                                printf("Error!! No 128data\n");
                            }
                            break;                           
                        }
                        fscanf(input_file,"%d,%d,%d,%d,%d,%d,%d,%d,%d,",&No_con_xi_1[l],&No_con_eta_1[l],&Order_xi[l],&Order_eta[l],&number5_a,&number5_b,&number5_c,&number5_d,&number5_e);
                        No_con_xi[l]=No_con_xi_1[l]+1;
                        No_con_eta[l]=No_con_eta_1[l]+1;
                        printf("%d\n",number1_a);
                        printf("No_con = (%d,%d)\n",No_con_xi[l],No_con_eta[l]);
                        printf("Order = (%d,%d)\n",Order_xi[l],Order_eta[l]);
                        No_knot_xi[l]=No_con_xi[l]+Order_xi[l]+1;
                        No_knot_eta[l]=No_con_eta[l]+Order_eta[l]+1;
                        printf("No_knot = (%d,%d)\n",No_knot_xi[l],No_knot_eta[l]);
                        all_number[l]=No_knot_xi[l]+No_knot_eta[l]+No_con_xi[l]*No_con_eta[l]*4+4;
                        //printf("%d%d%d%d%d\n",number5_a,number5_b,number5_c,number5_d,number5_e);
                        j=0;
                        do{
                            fscanf(input_file,"%lf%c%d,",&number,&D,&power);
                            //printf("%lf %d %d\n",number,D,power);
                            D_number[l][j]=number*pow(10,power);
                            //printf("D_number[%d] = %lf\n",j,D_number[l][j]);
                            if(D == 68)
                            {
                                j++;
                            }
                            else
                            {
                                j=j;
                            }
                        }
                        while(j < all_number[l]);
            
                        k=0;
                        for(i=0;i<No_knot_xi[l];i++){
                            knot_xi[l][i]=D_number[l][k];
                            printf("knot_xi[%d]=%lf\n",i,knot_xi[l][i]);
                            k++;
                        }
                        printf("\n");
                        for(ii=0;ii<No_knot_eta[l];ii++){
                            knot_eta[l][ii]=D_number[l][k];
                            printf("knot_eta[%d]=%lf\n",ii,knot_eta[l][ii]);
                            k++;
                        }
                        printf("\n");
                        for(iii=0;iii<No_con_xi[l]*No_con_eta[l];iii++){
                            weight[l][iii]=D_number[l][k];
                            printf("weight[%d]=%lf\n",iii,weight[l][iii]);
                            k++;
                        }
                        printf("\n");
                        for(iiii=0;iiii<No_con_xi[l]*No_con_eta[l];iiii++){
                            for(jjjj=0;jjjj<3;jjjj++){
                                con_coord[l][iiii][jjjj]=D_number[l][k];
                                printf("con_coord[%d][%d]=%lf\n",iiii,jjjj,con_coord[l][iiii][jjjj]);
                                k++;
                            }
                        }
                        ll++;
                    }
                    else{
                        //printf("No char_P!\n");
                    }
                }
            }
            printf("\n");
            break;
        }
    }
    fclose(input_file);
}

void make_input_iga(char output_file_name[])
{
    FILE *output_file;
    output_file=fopen(output_file_name,"w");
    assert(output_file != NULL);

    int i,j,k,l;

    fprintf(output_file,"E  nu\n\n");
    fprintf(output_file,"No_patch\n\n");
    fprintf(output_file,"Total_Control_Point\n\n");
    for(l=0;l<patch_no;l++){
        fprintf(output_file,"%d\t%d\n",Order_xi[l],Order_eta[l]);
    }
    fprintf(output_file,"\n");
    for(l=0;l<patch_no;l++){
        fprintf(output_file,"%d\t%d\n",No_knot_xi[l],No_knot_eta[l]);
    }
    fprintf(output_file,"\n");
    for(l=0;l<patch_no;l++){
        fprintf(output_file,"%d\t%d\n",No_con_xi[l],No_con_eta[l]);
    }
    fprintf(output_file,"\n");
    fprintf(output_file,"Patch_connectivity\n\n");

    fprintf(output_file,"Constraint\tLoad\tDistributedForce\n\n");

    for(l=0;l<patch_no;l++){
        for(i=0;i<No_knot_xi[l];i++){
            fprintf(output_file,"%.10lf ",knot_xi[l][i]);
        }
        fprintf(output_file,"\n");
        for(j=0;j<No_knot_eta[l];j++){
            fprintf(output_file,"%.10lf ",knot_eta[l][j]);
        }
        fprintf(output_file,"\n");
    }
    fprintf(output_file,"\n");

    for(l=0;l<patch_no;l++){
        for(k=0;k<No_con_xi[l]*No_con_eta[l];k++){
            fprintf(output_file,"%d\t%.20lf %.20lf %.20lf\n",k,con_coord[l][k][0],con_coord[l][k][1],weight[l][k]); //ここいじった
        }
        fprintf(output_file,"\n");
    }

    fprintf(output_file,"Constraint\n\nLoad\n\nDistributedForce\n");

    fclose(output_file);
}

int main(int argc, char *argv[])
{
    char input_file_name[80];
    char output_file_name[80];

    strcpy(input_file_name,argv[1]);
    strcpy(output_file_name,argv[2]);
    printf("--------------------------------------\n");
    printf("input_file_name [.igs] = %s\n",input_file_name);
    printf("output_file_name [.txt] = %s\n",output_file_name);
    printf("--------------------------------------\n");

    read_igs(input_file_name);
    make_input_iga(output_file_name);

    return(0);
}