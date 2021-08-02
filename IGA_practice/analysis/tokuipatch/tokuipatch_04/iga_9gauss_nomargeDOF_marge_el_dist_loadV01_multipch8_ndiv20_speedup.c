#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define ERROR					-999
//#define ERROR					0


#define	MAX_NO_CCpoint_ON_ELEMENT		16						//分割節点数
#define DIMENSION				2						//次元数
#define MAX_KIEL_SIZE				MAX_NO_CCpoint_ON_ELEMENT*DIMENSION	//要素分割マトリックスの大きさ
#define Ng						3						//Gauss-Legendreの足す回数
#define POW_Ng					Ng * Ng			//NgのDIMENSION乗の計算
#define	D_MATRIX_SIZE			3						//応力歪マトリックスの大きさ（2次元:3 3次元:6）

#define K_DIVISION_LENGE		10					//全体剛性マトリックスのcol&ptrを制作時に分ける節点数
#define EPS						0.0000000001				//連立1次方程式の残差
#define	N_STRAIN				4
#define N_STRESS				4
//各種最大配置可能数
#define MAX_N_KNOT 				1000
#define MAX_N_ELEMENT			110000
#define MAX_N_NODE				110000
#define	MAX_N_LOAD				100000
#define MAX_N_CONSTRAINT		100000
#define MAX_K_WHOLE_SIZE		MAX_N_NODE*DIMENSION
#define MAX_NON_ZERO			10000000
#define MAX_N_PATCH 	100
#define MAX_N_Controlpoint_in_Patch 10000

#define MAX_N_DISTRIBUTE_FORCE	100
#define DISTRIBUTE_FORCE_Ng		3
#define element_ndiv 20
//void Force_Dis( int Total_DistributeForce, int DistributeForce[MAX_N_DISTRIBUTE_FORCE][3], double Val_DistributeForce[MAX_N_DISTRIBUTE_FORCE],int *Total_Load,int Load_Node_Dir[MAX_N_LOAD][2],double Value_of_Load[MAX_N_LOAD],int Total_Control_Point, int El_No, int *Total_Element );
int Make_K_EL(int El_No, double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], double K_EL[MAX_KIEL_SIZE][MAX_KIEL_SIZE], double E, double nu ,int DM ,  int Total_Element, int Total_Control_Point);
void Get_InputData(	double *E,double *nu,int *Total_Element, int *Total_Control_Point,
			int *Total_Load,int *No_Patch,int Load_Node_Dir[MAX_N_LOAD][2],double Value_of_Load[MAX_N_LOAD],
			int *Total_Constraint, int Constraint_Node_Dir[MAX_N_CONSTRAINT][2],double Value_of_Constraint[MAX_N_CONSTRAINT],
			int *Total_DistributeForce, int argc, char *argv[]);
//全体剛性マトリックス

int Make_Index_Dof( int Total_Control_Point, int Total_Constraint, int Constraint_Node_Dir[MAX_N_CONSTRAINT][2] );
void Make_K_Whole_Ptr_Col(int Total_Element, int Total_Control_Point, int K_Whole_Size );
void Make_K_Whole_Val( double E, double nu, int Total_Element, int K_Whole_Size, int DM, int Total_Control_Point/*,int real_element[MAX_N_ELEMENT]*/);
//連立1次方程式
void Make_F_Vec( int Total_Load,int Load_Node_Dir[MAX_N_LOAD][2],double Value_of_Load[MAX_N_LOAD],int K_Whole_Size );
void Make_F_Vec_disp_const(int Total_Constraint,int Constraint_Node_Dir[MAX_N_CONSTRAINT][2],double Value_of_Constraint[MAX_N_CONSTRAINT],int Total_Element,double  E,double nu, int DM, int Total_Control_Point);
void mat_vec_crs(double vec_result[], double vec[], const int ndof);
double inner_product(int ndof, double vec1[], double vec2[]);
int check_conv_CG(int ndof, double alphak, double pp[], double eps, int itr);
void Diag_Scaling_CG_pre(int ndof, int flag_operation);
void CG_Solver(int ndof, int max_itr, double eps, int flag_ini_val);
//各種値
void Make_Strain(double E, double nu, int Total_Element, int El_No,int Total_Control_Point);
void Make_Stress_2D(double E, double nu, int Total_Element,int DM);
void Make_ReactionForce( int Total_Element,int Total_Control_Point,int El_No);
void Make_Parameter_z(int Total_Element, double E,double nu, int DM);
//分布荷重
void Force_dis(int Distriction_Force[DIMENSION][3], double Val_Distribute_Force[DIMENSION], double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], double Fe[DIMENSION] );
//void Make_Output( int Total_Control_Point, int Total_Element );
//NURBSの計算
void element_coordinate(int Total_Element,int Total_Control_Point);
void calculate_Controlpoint_using_NURBS(double element[DIMENSION],int Total_Element,int Total_Control_Point,double E, double nu,int DM);
void calculate_extendmesh_using_NURBS(double element_emsh[DIMENSION],int Total_Element,int Total_Control_Point);
void Gausspoint_coordinate(int Total_Element,int Total_Control_Point);
//void Make_Stress_2D_refine(double E, double nu, int Total_Element,int DM,int Total_Control_Point) ;
/*///J積分
int Make_B_x_Matrix_Quad_4(double B_x[DIMENSION][KIEL_SIZE], double Local_coord[DIMENSION], double X[No_Control_point_ON_ELEMENT][DIMENSION], double *J );
void Make_Strain_x_Quad_4(double E, double nu, int Total_Element);
void Make_EMT(double E, double nu, int Total_Element);*/

/* Distributed Load */

int  SerchForElement(int iPatch, int Total_Element, int iX, int iY);

void Setting_Dist_Load_2D(int Total_Control_Point, int iPatch, int Total_Element, int iCoord, double val_Coord,
         double Range_Coord[2], int type_load, double Coeff_Dist_Load[3]);

void Add_Equivalent_Nodal_Forec_to_F_Vec(int Total_Control_Point);


//static int DIMENSION;
static int KIEL_SIZE;//要素分割マトリックスの大きさ
//int No_Control_point_ON_ELEMENT=1;


static int Controlpoint_of_Element_nomerge[MAX_N_ELEMENT][MAX_NO_CCpoint_ON_ELEMENT];
static int Controlpoint_of_Element[MAX_N_ELEMENT][MAX_NO_CCpoint_ON_ELEMENT];
static double Node_Coordinate[MAX_N_NODE][DIMENSION+1];
static double Equivalent_Nodal_Force[MAX_N_NODE][DIMENSION]; // Equivalent nodal forces arising from the distributed load
static int K_Whole_Ptr[MAX_K_WHOLE_SIZE+1],  K_Whole_Col[MAX_NON_ZERO];
static double K_Whole_Val[MAX_NON_ZERO];
static int Index_Dof[MAX_K_WHOLE_SIZE];
static int INC[MAX_N_PATCH][MAX_N_ELEMENT][DIMENSION];
static int Adress_Controlpoint[MAX_N_PATCH][1000][1000]; //INCの配列をいじったものAdress_Controlpoint[ξ][η]；コントールポイント番号
static int Order[MAX_N_PATCH][DIMENSION];
static int No_knot[MAX_N_PATCH][DIMENSION];
static int No_Control_point[MAX_N_PATCH][DIMENSION];
static double element_coordinate_Nopoint[MAX_N_ELEMENT][DIMENSION];
static double Gausspoint_coordinates[MAX_N_ELEMENT][POW_Ng][DIMENSION];
//static int same_point[100];
static int same_point_in_Element[MAX_N_NODE];
static int Patch_controlpoint[MAX_N_PATCH][MAX_N_Controlpoint_in_Patch];//バッチとコントロールポイント番号の要素コネクティビティ
static int Element_patch[MAX_N_ELEMENT];//要素がどのパッチに属しているか示す配列(要素番号は1つのモデルで通し番号)
static int No_Controlpoint_in_patch[MAX_N_PATCH];
static int No_Control_point_ON_ELEMENT[10000];

static int Node_To_Node[K_DIVISION_LENGE][10000], Total_Control_Point_To_Node[K_DIVISION_LENGE];//ある節点に関係する節点番号s
//static int col_N[10][1000];

static double sol_vec[MAX_K_WHOLE_SIZE];
static double rhs_vec[MAX_K_WHOLE_SIZE];
static double diag_scaling[MAX_K_WHOLE_SIZE];

static double Shape[DIMENSION][MAX_N_NODE][10];
static double shape_func[MAX_N_NODE];
static double dShape_func1[MAX_N_NODE];
static double dShape_func2[MAX_N_NODE];
static double dShape[DIMENSION][MAX_N_NODE];
static double Position_Knots[MAX_N_PATCH][DIMENSION][MAX_N_KNOT];
static double Position_Data_param[DIMENSION];

static double Displacement[MAX_K_WHOLE_SIZE];
static double Strain[MAX_N_ELEMENT][POW_Ng][N_STRAIN];
static double Stress[MAX_N_ELEMENT][POW_Ng][N_STRESS];
static double ReactionForce[MAX_K_WHOLE_SIZE];

static double difference[MAX_N_PATCH][MAX_N_KNOT][DIMENSION];    /*隣り合うノットベクトルの差*/
static int ENC[MAX_N_ELEMENT][DIMENSION];   /*ENC[全ての要素][0,1]=x,y方向の何番目の要素か*/
static int real_Total_Element;	/*ゼロエレメントの要素数*/
static int real_element[MAX_N_ELEMENT];	 /*ゼロエレメントではない要素の番号*/
static int Total_element_all_ID[MAX_N_ELEMENT];    /*ゼロエレメントではない要素＝１、ゼロエレメント＝０*/
static int line_No_Total_element[MAX_N_PATCH][DIMENSION];   /*ゼロエレメントを含むすべての要素列の数*/
static int line_No_real_element[MAX_N_PATCH][DIMENSION];   /*ゼロエレメントではない要素列の数*/
static int real_element_line[MAX_N_PATCH][MAX_N_ELEMENT][DIMENSION];   /*ゼロエレメントではない要素列*/

static int No_points_for_colored_points;   /*zarusobaで点に色付ける時の全ての点の数*/
//static double data_result_shape_x[100000];
//static double data_result_shape_y[100000];
//static double data_result_disp_x[100000];
//static double data_result_disp_y[100000];

static int No_points_for_new_zarusoba;   /*zarusobaで点に色付ける時の全ての点の数*/
static double data_result_shape_x_for_new_zarusoba[300000];
static double data_result_shape_y_for_new_zarusoba[300000];
static double data_result_disp_x_for_new_zarusoba[300000];
static double data_result_disp_y_for_new_zarusoba[300000];



//static double Strain_x[MAX_N_ELEMENT][POW_Ng][N_STRAIN];


FILE *fp;

int main(int argc, char *argv[]){
	clock_t start,end;
	//printf("aaaaaaa\n");
	int i,j, k,l;
	int re;
  //int p;
	int q,r;
	int DM=1;
	int Total_Element;
	int Total_Control_Point;
	int No_Patch=0;
	int Total_net=0;
	static int Total_Load=0,			Load_Node_Dir[MAX_N_LOAD][2];		static double Value_of_Load[MAX_N_LOAD];
	static int Total_Constraint=0,	Constraint_Node_Dir[MAX_N_CONSTRAINT][2];	static double Value_of_Constraint[MAX_N_CONSTRAINT];
	static int Total_DistributeForce=0;
	int K_Whole_Size=0;
	int El_No=0;
	static double element[DIMENSION];
	static double element_emsh[DIMENSION];
	static double E,nu;
	static int max_itr;

	if (argc != 2) {
		printf("Argument is missing\n");
	}


	//printf("bbbbb\n");
	start = clock();
	Get_InputData( &E,&nu, &Total_Element, &Total_Control_Point, &Total_Load,&No_Patch
		       ,Load_Node_Dir, Value_of_Load, &Total_Constraint, Constraint_Node_Dir, Value_of_Constraint, &Total_DistributeForce, argc, argv);
	printf("Finish Get_InputData\n");
	printf("Total Element=%d Node=%d Constraint=%d Load=%d\n",Total_Element ,Total_Control_Point, Total_Constraint, Total_Load);
	printf("E;%le nu;%le\n", E,nu);

	/////////////////////全体剛性マトリックスの制作////////////////////////////
	K_Whole_Size = Make_Index_Dof( Total_Control_Point, Total_Constraint, Constraint_Node_Dir );
	//printf("K_Whole_Size=%d\n",K_Whole_Size);
	Make_K_Whole_Ptr_Col( Total_Element, Total_Control_Point, K_Whole_Size );
	Make_K_Whole_Val( E, nu, Total_Element, K_Whole_Size, DM, Total_Control_Point/*,real_element[MAX_N_ELEMENT]*/);
	printf("Finish Make_K_Whole\n");
	printf("check\n");
	///////////////連立一次方程式/////////////////////////////////////////
//printf("DistributeForce=%d\n", DistributeForce);
	//printf("Val_DistributeForce=%d\n", Val_DistributeForce);
	//printf("Load_Node_Dir=%s\n", );
	/*		Force_Dis(Total_DistributeForce, DistributeForce, Val_DistributeForce, &Total_Load, Load_Node_Dir, Value_of_Load, Total_Control_Point,El_No,&Total_Element);*/

	for(i=0; i<Total_Load; i++)
		printf("%11.10e\n",Value_of_Load[i]);
		printf("pp");
	max_itr = K_Whole_Size;
	printf("K_Whole_Size:%d\n",K_Whole_Size);
	Make_F_Vec( Total_Load, Load_Node_Dir, Value_of_Load, K_Whole_Size );
	Make_F_Vec_disp_const( Total_Constraint, Constraint_Node_Dir, Value_of_Constraint, Total_Element, E, nu, DM, Total_Control_Point);
	Add_Equivalent_Nodal_Forec_to_F_Vec(Total_Control_Point);

	Diag_Scaling_CG_pre( K_Whole_Size, 0);
	CG_Solver( K_Whole_Size, max_itr, EPS, 0);
	Diag_Scaling_CG_pre( K_Whole_Size,1);
	printf("Finish CG_Solver\n");
	/////////////変位と歪と応力//////////////////////////////////////
	for( i = 0; i < Total_Constraint; i++ )
		Displacement[ Constraint_Node_Dir[i][0]*DIMENSION + Constraint_Node_Dir[i][1] ] = Value_of_Constraint[i];
	for( i = 0; i < Total_Control_Point; i++ ){for( j = 0; j < DIMENSION; j++ ){
		int index = Index_Dof[i*DIMENSION+j];
		if( index >= 0 ) Displacement[i*DIMENSION+j] = sol_vec[index];
//		printf("Displacement[%d] = %le\n",i*DIMENSION+j,Displacement[i*DIMENSION+j]);
	}}
	printf("Finish Make_Displacement\n");
	end = clock();
	printf("%.2f秒かかりました\n",(double)(end-start)/CLOCKS_PER_SEC);
	Make_Strain( E, nu, Total_Element, El_No,Total_Control_Point);
	printf("Finish Make_Strain\n");
	Make_Stress_2D( E, nu, Total_Element,DM);
	printf("Finish Make_Stress\n");
	Make_ReactionForce( Total_Element, Total_Control_Point, El_No);
	printf("Finish Make_ReactionForce\n");
	puts("sol_vec");
	Make_Parameter_z( Total_Element, E, nu, DM);
	printf("Finish Make_Parameter_z\n");
	//Make_Output( Total_Control_Point, Total_Element );
	//printf("Finish Make_Output\n" );

	//Make_Strain_x_Quad_4( E, nu, Total_Element);
	//printf("Finish Make_Strain_x_Quad_4\n");
	//Make_EMT(E, nu, Total_Element);
	//printf("Finish Make_EMT\n");
/*

	for(i = 0; i <K_Whole_Size; i++ ){
		printf("%le ",sol_vec[i]);
		printf("\n");
	}
	puts("DIMENSION");
	for(j = 0; j < Total_Control_Point; j++ ){
		for(i = 0; i < DIMENSION; i++ )
			printf("%le ",Displacement[j*DIMENSION+i]);
		printf("\n");
	}
	printf("\n");
	puts("Strain");
	for( k = 0; k < POW_Ng; k++){
		for( i = 0; i < Total_Element; i++ ){
			for( j = 0; j < N_STRAIN; j++ )
				printf("%le ",Strain[i][k][j]);
			printf("\n");
		}
	}
	printf("\n");
	puts("Stress");
	for( k = 0; k < POW_Ng; k++){
		for( i = 0; i < Total_Element; i++ ){
			for( j = 0; j < N_STRESS; j++ )
				printf("%le ",Stress[i][k][j]);
			printf("\n");
		}
	}
	printf("\n");
	puts("ReactionForce");
	for(j = 0; j < Total_Control_Point; j++ ){
		for(i = 0; i < DIMENSION; i++ )
			printf("%le ",ReactionForce[j*DIMENSION+i]);
		printf("\n");
	}

*/




	fp = fopen("checkAns/checkAns.txt", "w");
	for(i = 0; i <K_Whole_Size; i++ ){
		fprintf(fp,"%le ",sol_vec[i]);
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n\n\nDisplacement\n");
	for(j = 0; j < Total_Control_Point; j++ ){
		fprintf(fp,"%d\t",j );
		for(i = 0; i < DIMENSION; i++ )
			fprintf(fp,"%.13e\t ",Displacement[j*DIMENSION+i]);
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n\n\nStrain\n");
	//for( i = 0; i < Total_Element; i++ ){
	for( re = 0; re < real_Total_Element; re++ ){
		i=real_element[re];
		for( k = 0; k < POW_Ng; k++){
			fprintf(fp,"%d\t%d\t",i,k);
			for( j = 0; j < N_STRAIN; j++ )
				fprintf(fp,"%.13e\t",Strain[i][k][j]);
			fprintf(fp,"\n");
		}
	}
/*	fprintf(fp,"\n\n\n\n\n\n\n\n\n\n\nStrain_x\n");
	for( k = 0; k < POW_Ng; k++){
		for( i = 0; i < Total_Element; i++ ){
			for( j = 0; j < DIMENSION; j++ )
				fprintf(fp,"%d	%d	%le\t ",k,i,Strain_x[i][k][j]);
			fprintf(fp,"\n");
		}
	}*/

	fprintf(fp,"\n\n\n\n\n\n\n\n\n\n\nStress\n");
	//for( i = 0; i < Total_Element; i++ ){
	for(re=0;re<real_Total_Element;re++){
		i=real_element[re];
		for( k = 0; k < POW_Ng; k++){
			fprintf(fp,"%d\t%d\t",i,k);
			for( j = 0; j < N_STRESS; j++ )
				fprintf(fp,"%.13e\t",Stress[i][k][j]);
			fprintf(fp,"\n");
		}
	}
	fprintf(fp,"\n\n\nReaction Force\n");
	for(j = 0; j < Total_Control_Point; j++ ){
		for(i = 0; i < DIMENSION; i++ )
			fprintf(fp,"%.13e\t ",ReactionForce[j*DIMENSION+i]);
		fprintf(fp,"\n");
	}
	fclose(fp);

	fp = fopen("Displacement.dat", "w");
	fprintf(fp,"label=Displacement\n");
	fprintf(fp, "num_items=%d\n", Total_Control_Point);
	fprintf( fp, "\n");
	for(j = 0; j < Total_Control_Point; j++ ){
			fprintf(fp,"%d:	%le %le ",j,Displacement[j*DIMENSION+0],Displacement[j*DIMENSION+1]);
		fprintf(fp,"\n");
	}
	fclose(fp);

	fp = fopen("Strain@IntegrationPoint.dat", "w");
	fprintf(fp,"label=Strain@IntegrationPoint\n");
	fprintf(fp, "num_items=%d\n", real_Total_Element);
	fprintf( fp, "\n");
	//for( i = 0; i < Total_Element; i++ ){
	for(re=0;re<real_Total_Element;re++){
		i=real_element[re];
		fprintf(fp,"%d:\t", i);
		for( k = 0; k < POW_Ng; k++){
				fprintf(fp,"%13e\t%13e\t%13e\t%13e\t",Strain[i][k][0],Strain[i][k][1],Strain[i][k][2],Strain[i][k][3]);
			}fprintf(fp,"\n");
		}
	fclose(fp);

	fp = fopen("Stress@IntegrationPoint.dat", "w");
	fprintf(fp,"label=Stress@IntegrationPoint\n");
	fprintf(fp, "num_items=%d\n", real_Total_Element);
	fprintf( fp, "\n");

	//for( i = 0; i < Total_Element; i++ ){
	for(re=0;re<real_Total_Element;re++){
		i=real_element[re];
		fprintf(fp,"%d\t:",i );
		for( k = 0; k < POW_Ng; k++){
			fprintf(fp,"%13e\t%13e\t%13e\t%13e\t",Stress[i][k][0],Stress[i][k][1],Stress[i][k][2],Stress[i][k][3]);
		}
		fprintf(fp,"\n");
	}

	fclose(fp);


	fp = fopen("ReactionForce.dat", "w");
	fprintf(fp,"label=ReactionForce\n");
	fprintf(fp, "num_items=%d\n", Total_Control_Point);
	fprintf( fp, "\n");
	for(j = 0; j < Total_Control_Point; j++ ){
			fprintf(fp,"%d:	%13e %13e ",j,ReactionForce[j*DIMENSION+0],ReactionForce[j*DIMENSION+1]);
		fprintf(fp,"\n");
	}
	fclose(fp);


	fp = fopen("mesh.msh", "w");
	fprintf(fp,"%d\n",real_Total_Element/*Total_Element*/);
	//for( i = 0; i < Total_Element; i ++ ){
	for(re=0;re<real_Total_Element;re++){
		i=real_element[re];
		fprintf(fp,"%d %d %d %d %d %d %d %d %d",Controlpoint_of_Element[i][8],Controlpoint_of_Element[i][6],Controlpoint_of_Element[i][0],Controlpoint_of_Element[i][2],Controlpoint_of_Element[i][7],Controlpoint_of_Element[i][3],Controlpoint_of_Element[i][1],Controlpoint_of_Element[i][5],Controlpoint_of_Element[i][4] );
		//fprintf(fp,"%d %d %d %d %d %d %d %d %d",Controlpoint_of_Element[i][0],Controlpoint_of_Element[i][2],Controlpoint_of_Element[i][8],Controlpoint_of_Element[i][6],Controlpoint_of_Element[i][1],Controlpoint_of_Element[i][5],Controlpoint_of_Element[i][7],Controlpoint_of_Element[i][3],Controlpoint_of_Element[i][4] );
		fprintf(fp, "\n" );
	}

	fprintf(fp, "%d\n", Total_Control_Point);
	for(i = 0; i < Total_Control_Point; i++ ){
		fprintf(fp, "%lf %lf\n", Node_Coordinate[i][0],Node_Coordinate[i][1]);
	}
	fprintf( fp, "\n");
	fclose(fp);


	fp=fopen("mesh_net/control_net.msh", "w");
	for (l = 0; l < No_Patch; l++) {
		Total_net += (No_Control_point[l][0]-1)*(No_Control_point[l][1]-1);
	}

	fprintf(fp, "%d\n",Total_net);
	for ( l = 0; l < No_Patch; l++) {
		for (j = 0; j <No_Control_point[l][1]-1; j++) {
					for (i = 0; i < No_Control_point[l][0]-1; i++) {
						fprintf(fp, "%d %d %d %d\n", Adress_Controlpoint[l][i][j],Adress_Controlpoint[l][i+1][j],Adress_Controlpoint[l][i+1][j+1],Adress_Controlpoint[l][i][j+1]);
					}
				}
	}

		fprintf(fp, "%d\n", Total_Control_Point);

		for(i = 0; i < Total_Control_Point; i++ ){
		fprintf(fp, "%lf %lf\n", Node_Coordinate[i][0],Node_Coordinate[i][1]);
		}
		fclose(fp);


	fp=fopen("mesh_net/element.msh","w");

	fprintf(fp, "%d\n",Total_Element);
 	element_coordinate(Total_Element,Total_Control_Point);
	for ( i= 0; i < Total_Element*9; i+=9) {
			fprintf(fp, "%d %d %d %d %d %d %d %d %d\n",same_point_in_Element[i],same_point_in_Element[i+1],same_point_in_Element[i+2],same_point_in_Element[i+3],same_point_in_Element[i+4],
			same_point_in_Element[i+5],same_point_in_Element[i+6],same_point_in_Element[i+7],same_point_in_Element[i+8]);
		}
		fprintf(fp, "%d\n", Total_Element*9);
		for (i = 0; i < Total_Element*9; i++) {
			fprintf(fp,"%lf %lf\n",element_coordinate_Nopoint[i][0],element_coordinate_Nopoint[i][1]);
			}
	fclose(fp);

/*(2019_10_10)
	fp = fopen("colored_point/NURBS_stress_x.dat","w");
	fprintf(fp,"label=Stress x \nnum_items=%d\n\n",(element_ndiv+1)*(element_ndiv+1)*real_Total_Element);
	fclose(fp);

	fp = fopen("colored_point/NURBS_stress_y.dat","w");
	fprintf(fp,"label=Stress y \nnum_items=%d\n\n",(element_ndiv+1)*(element_ndiv+1)*real_Total_Element);
	fclose(fp);

	fp = fopen("colored_point/NURBS_stress_xy.dat","w");
	fprintf(fp,"label=Stress xy \nnum_items=%d\n\n",(element_ndiv+1)*(element_ndiv+1)*real_Total_Element);
	fclose(fp);

	fp = fopen("colored_point/NURBS_stress_z.dat","w");
	fprintf(fp,"label=Stress z \nnum_items=%d\n\n",(element_ndiv+1)*(element_ndiv+1)*real_Total_Element);
	fclose(fp);

	fp = fopen("colored_point/NURBS_points.txt","w");
	fprintf(fp,"%d\n",(element_ndiv+1)*(element_ndiv+1)*real_Total_Element);
	for(p=0;p<No_points_for_colored_points;p++){
			fprintf(fp,"%-.13lf  %-.13lf\n",data_result_shape_x[p],data_result_shape_y[p]);
	}

	fclose(fp);


	fp = fopen("colored_point/NURBS_disp_x.dat","w");

	fprintf(fp,"label=Displacement x\nnum_items=%d\n\n",(element_ndiv+1)*(element_ndiv+1)*real_Total_Element);
	  for(p=0;p<No_points_for_colored_points;p++){
			fprintf(fp,"%d:%-.13le\n",p,data_result_disp_x[p]);
	}

	fclose(fp);


	fp = fopen("colored_point/NURBS_disp_y.dat","w");

	fprintf(fp,"label=Displacement y\nnum_items=%d\n\n",(element_ndiv+1)*(element_ndiv+1)*real_Total_Element);
	for(p=0;p<No_points_for_colored_points;p++){
			fprintf(fp,"%d:%-.13le\n",p,data_result_disp_y[p]);
	}

	fclose(fp);(2019_10_10)*/

	//fp = fopen("colored_point/NURBS_disp_radius.dat","w");

	//fprintf(fp,"label=Displacement\nnum_items=%d\n\n",(element_ndiv+1)*(element_ndiv+1)*real_Total_Element);
	/*for(p=0;p<No_points_for_colored_points;p++){
			fprintf(fp,"%d:%-.13le\n",p,pow((data_result_disp_x[p]*data_result_disp_x[p]+data_result_disp_y[p]*data_result_disp_y[p]),0.5));
	}*/

	//fclose(fp);


	//fp = fopen("NURBS/NURBS_disp.dat", "w");

		calculate_Controlpoint_using_NURBS(element,Total_Element,Total_Control_Point,E,nu,DM);

		//fclose(fp);

		/*fp = fopen("NURBS/control_point.dat", "w");
		for(j = 0; j < Total_Control_Point; j++ ){
				fprintf(fp,"%d:  %le  %le	 %le %le ",j,Node_Coordinate[j][0],Node_Coordinate[j][1],Displacement[j*DIMENSION+0],Displacement[j*DIMENSION+1]);
			fprintf(fp,"\n");
		}
		fclose(fp);*/

		fp = fopen("new_zarusoba/control_point.dat", "w");
		fprintf(fp, "%d\n", Total_Control_Point);
		for(j = 0; j < Total_Control_Point; j++ ){
				fprintf(fp,"%le  %le\n",Node_Coordinate[j][0],Node_Coordinate[j][1]);
		}
		fclose(fp);

		calculate_extendmesh_using_NURBS(element_emsh,Total_Element,Total_Control_Point);

		fp = fopen("new_zarusoba/extended_mesh.emsh", "w");
		fprintf(fp,"%d\n",real_Total_Element*10);

		q=0;r=0;
		for (i = Total_Control_Point; i < (Total_Control_Point+No_points_for_new_zarusoba-2); i=i+2) {
			if (q != 10+11*r) {
			fprintf(fp,"%d %d %d %d\n", i,i+1,i+3,i+2);
			//printf("q:%d\n",q );
		}
		if (q==10+11*r) {
			r++;
		}
			//printf("r:%d\n",r);
			q++;
		}

		for (l = 0; l < No_Patch; l++) {
		Total_net += (No_Control_point[l][0]-1)*(No_Control_point[l][1]-1);
	}

	fprintf(fp, "%d\n",Total_net*2);
	for ( l = 0; l < No_Patch; l++) {
		for (j = 0; j <No_Control_point[l][1]-1; j++) {
					for (i = 0; i < No_Control_point[l][0]-1; i++) {
						fprintf(fp, "%d %d\n",Adress_Controlpoint[l][i][j],Adress_Controlpoint[l][i+1][j]);
						fprintf(fp, "%d %d\n",Adress_Controlpoint[l][i+1][j],Adress_Controlpoint[l][i+1][j+1]);
						fprintf(fp, "%d %d\n",Adress_Controlpoint[l][i+1][j+1],Adress_Controlpoint[l][i][j+1]);
						fprintf(fp, "%d %d\n",Adress_Controlpoint[l][i][j+1],Adress_Controlpoint[l][i][j]);
					}
				}
	}

		fprintf(fp, "%d\n", Total_Control_Point+No_points_for_new_zarusoba);

		for(i = 0; i < Total_Control_Point; i++ ){
		fprintf(fp, "%lf %lf\n", Node_Coordinate[i][0],Node_Coordinate[i][1]);
		}

		for ( i = Total_Control_Point; i < Total_Control_Point+No_points_for_new_zarusoba; i++) {
			fprintf(fp,"%lf %lf\n",data_result_shape_x_for_new_zarusoba[i],data_result_shape_y_for_new_zarusoba[i]);
		}

		fclose(fp);

		fp = fopen("new_zarusoba/radial_displacement.dat","w");
    fprintf(fp,"label=Radial Displacement\nnum_items=%d\n\n",Total_Control_Point+No_points_for_new_zarusoba);
		for (i = 0; i < Total_Control_Point; i++) {
			fprintf(fp,"%d:0.0\n",i);
		}
		for (i = Total_Control_Point; i < (Total_Control_Point+No_points_for_new_zarusoba); i++)
				fprintf(fp,"%d:%.16e\n",i,pow((data_result_disp_x_for_new_zarusoba[i]*data_result_disp_x_for_new_zarusoba[i]+data_result_disp_y_for_new_zarusoba[i]*data_result_disp_y_for_new_zarusoba[i]),0.5));
    fclose(fp);

		fp = fopen("new_zarusoba/displacement_y.dat","w");
		fprintf(fp,"label=Displacement y\nnum_items=%d\n\n",Total_Control_Point+No_points_for_new_zarusoba);
		for (i = 0; i < Total_Control_Point; i++) {
			fprintf(fp,"%d:%.16e\n",i,Displacement[i*DIMENSION+1]);
		}
		for (i = Total_Control_Point; i < (Total_Control_Point+No_points_for_new_zarusoba); i++)
				fprintf(fp,"%d:%.16e\n",i,data_result_disp_y_for_new_zarusoba[i]);
		fclose(fp);

		fp = fopen("new_zarusoba/displacement_x.dat","w");
		fprintf(fp,"label=Displacement x\nnum_items=%d\n\n",Total_Control_Point+No_points_for_new_zarusoba);
		for (i = 0; i < Total_Control_Point; i++) {
			fprintf(fp,"%d:%.16e\n",i,Displacement[i*DIMENSION+0]);
		}
		for (i = Total_Control_Point; i < (Total_Control_Point+No_points_for_new_zarusoba); i++)
				fprintf(fp,"%d:%.16e\n",i,data_result_disp_x_for_new_zarusoba[i]);
		fclose(fp);
	/*printf("Finish Make_Displacement\n");
	Make_Strain_Quad_4( E, nu, Total_Element, El_No,Total_Control_Point);
	printf("Finish Make_Strain\n");
	Make_Stress_2D( E, nu, Total_Element, DM);
	printf("Finish Make_Stress\n");
	Make_ReactionForce_Quad_4( Total_Element, Total_Control_Point, El_No );
	printf("Finish Make_ReactionForce\n");
	puts("sol_vec");
	Make_Parameter_z( Total_Element, E, nu, DM);
	printf("Finish Make_Parameter_z\n");
	Make_Output( Total_Control_Point, Total_Element );
	printf("Finish Make_Output\n" );*/






		fp = fopen("Gauss_stress/Gausspoint_coordinates.dat","w");
		fprintf(fp, "%d\n",Total_Element*9);
	 	Gausspoint_coordinate(Total_Element,Total_Control_Point);
			for (i = 0; i < Total_Element; i++) {
				for ( j = 0; j < POW_Ng; j++) {
					fprintf(fp,"%.16e %.16e\n",Gausspoint_coordinates[i][j][0],Gausspoint_coordinates[i][j][1]);
				}
				}
		fclose(fp);

		fp = fopen("Gauss_stress/Gauss_stress_x.dat","w");
		k=0;
    fprintf(fp,"label=Stress x\nnum_items=%d\n\n",Total_Element*9);
		for (i = 0; i < Total_Element; i++) {
			for ( j = 0; j < POW_Ng; j++) {
				fprintf(fp,"%d:%.16e\n",k,Stress[i][j][0]);
				k++;
			}
			}

    fclose(fp);

		fp = fopen("Gauss_stress/Gauss_stress_y.dat","w");
		k=0;
    fprintf(fp,"label=Stress y\nnum_items=%d\n\n",Total_Element*9);
		for (i = 0; i < Total_Element; i++) {
			for ( j = 0; j < POW_Ng; j++) {
				fprintf(fp,"%d:%.16e\n",k,Stress[i][j][1]);
				k++;
			}
			}

    fclose(fp);

		fp = fopen("Gauss_stress/Gauss_stress_xy.dat","w");

    fprintf(fp,"label=Stress xy\nnum_items=%d\n\n",Total_Element*9);
		k=0;
		for (i = 0; i < Total_Element; i++) {
			for ( j = 0; j < POW_Ng; j++) {
				fprintf(fp,"%d:%.16e\n",k,Stress[i][j][2]);
				k++;
			}
			}

    fclose(fp);



return 0;



}




///////////////////////////////////////////////////////
//////////////全体剛性マトリックスの制作///////////////
///////////////////////////////////////////////////////

//ファイルからデータをもらう
void Get_InputData(	double *E,double *nu,int *Total_Element, int *Total_Control_Point,
				   int *Total_Load,int *No_Patch,int Load_Node_Dir[MAX_N_LOAD][2],double Value_of_Load[MAX_N_LOAD],
				   int *Total_Constraint,int Constraint_Node_Dir[MAX_N_CONSTRAINT][2],double Value_of_Constraint[MAX_N_CONSTRAINT],
			int *Total_DistributeForce, int argc,char *argv[]){
  int i,j,k,l,iii;
	int n,p,q,h,x,y;
	char s[256];
	int ii,jj,kk,kkk;
	int e,b,B;
	int iiloc,jjloc,kkloc;
	int r=0;
	/* for the distributed loads*/


	if ((fp = fopen(argv[1], "r")) == NULL)	printf("file open error!!\n");
	//次元数
	/*fscanf(fp,"%d",&DIMENSION );
	printf("DIMENSION=%d\n",DIMENSION);
	fgets(s, 256, fp);*/
	//材料定数
		fscanf(fp, "%le %le", &*E, &*nu);
		fgets(s, 256, fp);
		printf("E:%le nu:%le\n", *E,*nu);
	//パッチ数
		fscanf(fp, "%d", &*No_Patch );
		fgets(s, 256, fp);
		printf("No_Patch:%d\n", *No_Patch);

		//コントロールポイント数
		fscanf(fp,"%d",&*Total_Control_Point);
		fgets(s, 256, fp);
		printf("Total_Control_Point:%d\n", *Total_Control_Point);
		//ξη方向の各次数
		for (l = 0; l < *No_Patch; l++) {
		for (j = 0; j < DIMENSION; j++) {
			fscanf(fp,"%d",&Order[l][j]);
			printf("Order[%d][%d]=%d\n",l,j,Order[l][j]);
		}
	}

		fgets(s, 256, fp);
		//ノット数
		for ( l = 0; l < *No_Patch; l++) {
			for (j = 0; j < DIMENSION; j++) {
				fscanf(fp,"%d",&No_knot[l][j]);
				printf("No_knot[%d][%d]=%d\n",l,j,No_knot[l][j]);
			}
		}

		fgets(s, 256, fp);
		//各パッチ各方向のコントロールポイント数
		for ( l = 0; l < *No_Patch; l++) {
			for (j = 0; j < DIMENSION; j++) {
				fscanf(fp,"%d",&No_Control_point[l][j]);
				printf("No_Control_point[%d][%d]:%d\n",l,j,No_Control_point[l][j] );
			}
		}

		fgets(s, 256, fp);

		for ( l = 0; l < *No_Patch; l++) {
			No_Controlpoint_in_patch[l]=1.0;
		}

		for (l = 0; l < *No_Patch; l++) {
			for ( j = 0; j < DIMENSION; j++) {
				No_Controlpoint_in_patch[l]*=No_Control_point[l][j];
			}
		}

		for ( l = 0; l < *No_Patch; l++) {
			for (j = 0; j < DIMENSION; j++) {
				if (No_knot[l][j]!=No_Control_point[l][j]+Order[l][j]+1) {
					printf("wrong relationship between the number of knot vector and the number of control_point \n");
					printf("in patch_No.%d direction:%d\n", l,j);
					exit(0);
				}
			}
		}


		for ( l = 0; l < *No_Patch; l++) {
			printf("No_Controlpoint_in_patch[%d]:%d\t",l,No_Controlpoint_in_patch[l]);
		}
		printf("\n");

		for (l = 0; l < *No_Patch; l++) {
			//printf("l;%d\n",l);
			for ( i = 0; i < No_Controlpoint_in_patch[l]; i++) {
				//printf("i:%d\n",i );
				fscanf(fp,"%d",&Patch_controlpoint[l][i]);
			}
		}


		/*for (l = 0; l < *No_Patch; l++) {
			//printf("l:%d\n",l);
			for ( i = 0; i < No_Controlpoint_in_patch[l]; i++) {
				//printf("i:%d\n", i);
				printf("Patch_controlpoint[%d][%d]:%d\n",l,i,Patch_controlpoint[l][i]);
			}printf("\n");
		}*/


		fscanf(fp,"%d %d %d",Total_Constraint,Total_Load,Total_DistributeForce);
		printf("Total_Constraint;%d\n",*Total_Constraint );
		printf("Total_Load;%d\n",*Total_Load);
		printf("Total_DistributedForce;%d\n",*Total_DistributeForce);
		fgets(s, 256, fp);

		//ノットベクトルの読み込み
		for ( l = 0; l < *No_Patch; l++) {
			for (j = 0; j < DIMENSION; j++) {
				for (k = 0; k < No_knot[l][j]; k++) {
					fscanf(fp,"%le",&Position_Knots[l][j][k]);
					//printf("%le\t",Position_Knots[l][j][k]);
				}//printf("\n");
			}
		}

		for ( l = 0; l < *No_Patch; l++) {
			No_Control_point_ON_ELEMENT[l]=1.0;
		}

		*Total_Element =0.0;

		for (l = 0; l < *No_Patch; l++) {
			if (DIMENSION ==2) {
				*Total_Element += (No_Control_point[l][0]-Order[l][0])*(No_Control_point[l][1]-Order[l][1]);
				No_Control_point_ON_ELEMENT[l] = (Order[l][0]+1)*(Order[l][1]+1);
				}
				else {
					*Total_Element += (No_Control_point[l][0]-Order[l][0])*(No_Control_point[l][1]-Order[l][1]*(No_Control_point[l][2]-Order[l][2]));
					No_Control_point_ON_ELEMENT[l] =(Order[l][0]+1) *(Order[l][1]+1)*(Order[l][2]+1);
				}
			}
			printf("Total_Element=%d\n",*Total_Element);

			/*for ( l = 0; l < *No_Patch; l++) {
				printf("No_Control_point_ON_ELEMENT[%d]=%d\n",l,No_Control_point_ON_ELEMENT[l]);
			}*/

	//節点座標
	for(i = 0; i < *Total_Control_Point; i++ ){
		fscanf(fp, "%d", &ii );
		for( j = 0; j < DIMENSION+1; j++ )
			fscanf(fp,"%le",&Node_Coordinate[ii][j] );//Node_Coordinate[i][2]:重み
	}
	/*for(i = 0; i < *Total_Control_Point; i++ )
	for( j = 0; j < DIMENSION+1; j++ )printf("Node_Coordinate[%d][%d]=%e\n",i,j,Node_Coordinate[i][j] );*/
	fgets(s, 256, fp);
	//拘束
	for( i = 0; i < *Total_Constraint; i++ )
		fscanf(fp, "%d %d %le", &Constraint_Node_Dir[i][0], &Constraint_Node_Dir[i][1],&Value_of_Constraint[i]);
		for( i = 0; i < *Total_Constraint; i++ )printf("Constraint_Node_Dir[%d][0]= %d Constraint_Node_Dir[%d][1]=%d Value_of_Constraint[%d]= %e \n",i, Constraint_Node_Dir[i][0], i,Constraint_Node_Dir[i][1],i,Value_of_Constraint[i]);
	fgets(s, 256, fp);
	//荷重
	for( i = 0; i < *Total_Load; i++ )
		fscanf(fp, "%d %d %le", &Load_Node_Dir[i][0], &Load_Node_Dir[i][1],&Value_of_Load[i]);
		for( i = 0; i < *Total_Load; i++ )printf("Load_Node_Dir[%d][0]= %d Load_Node_Dir[%d][1]= %d Value_of_Load[%d]= %e\n",i,Load_Node_Dir[i][0],i,Load_Node_Dir[i][1],i,Value_of_Load[i] );

		int iPatch, iCoord, type_load;
		double Range_Coord[2], val_Coord, Coeff_Dist_Load[3];
		  int iPatch_array[MAX_N_DISTRIBUTE_FORCE], iCoord_array[MAX_N_DISTRIBUTE_FORCE], type_load_array[MAX_N_DISTRIBUTE_FORCE];
		double val_Coord_array[MAX_N_DISTRIBUTE_FORCE],Range_Coord_array[MAX_N_DISTRIBUTE_FORCE][2], Coeff_Dist_Load_array[MAX_N_DISTRIBUTE_FORCE][3];

	fgets(s, 256, fp);
		for(i=0; i < *Total_DistributeForce; i++){
		  fscanf(fp, "%d %d %d %lf %lf %lf %lf %lf %lf", &type_load, &iPatch, &iCoord, &val_Coord, &Range_Coord[0], &Range_Coord[1], &Coeff_Dist_Load[0], &Coeff_Dist_Load[1], &Coeff_Dist_Load[2]);
		  printf("Distibuted load nober: %d\n",i);
		  printf("type_load: %d  iPatch: %d iCoord: %d  val_Coord: %15.7e  Range_Coord: %15.7e  %15.7e\n Coef_Dist_Load: %15.7e %15.7e %15.7e\n",
			 type_load, iPatch, iCoord,
			 val_Coord,Range_Coord[0],Range_Coord[1],Coeff_Dist_Load[0],Coeff_Dist_Load[1],Coeff_Dist_Load[2]);
		/*
 		type_load: Direction of distributed load: 0-x direction, 1-y direction, 2-normal to the segemet/surface
		iPatch: Patch number to which the distributed load is assigned., 0, 1, ...
		iCoord: 0: Distributed load is applied to line along Xi axis.
                        1: Distributed load is applied to line along Eta axis
		val_Coord: その時のもう片方の座標
		Range_Coord[0]: Local coordinate value at which the distributed load starts.
		Range_Coord[1]: Local coordinate value at which the distributed load ends.
		Coeff_Dist_Load[0], &Coeff_Dist_Load[1], &Coeff_Dist_Load[2]: The coefficients of distributed load value:
			Coeff_Dist_Load[0]*Xi + Coeff_Dist_Load[1]*Xi + Coeff_Dist_Load[2]*Xi^2
		or
			Coeff_Dist_Load[0]*Xi + Coeff_Dist_Load[1]*Eta + Coeff_Dist_Load[2]*Eta^2
  		*/
		  type_load_array[i] = type_load;
		  iPatch_array[i] = iPatch;
		  iCoord_array[i] = iCoord;
		  val_Coord_array[i] = val_Coord;
		  Range_Coord_array[i][0] = Range_Coord[0];
		  Range_Coord_array[i][1] = Range_Coord[1];
		  Coeff_Dist_Load_array[i][0] =  Coeff_Dist_Load[0];
		  Coeff_Dist_Load_array[i][1] =  Coeff_Dist_Load[1];
		  Coeff_Dist_Load_array[i][2] =  Coeff_Dist_Load[2];
		}


		  /* Setting_Dist_Load_2D(&*Total_Control_Point,&iPatch, &*Total_Element, iCoord, val_Coord,
		     Range_Coord[2],  norm_dir, type_load, Coeff_Dist_Load); */
//自由度共有の計算
//(同じ座標を計算して要素コネクティビィティのコントロールポイント番号を入れ替える)
		/*for (ii = 0; ii < *Total_Control_Point; ii++) {
			same_point[ii]=ii;
		}

		for (ii = 0 ; ii < *Total_Control_Point; ii++) {
			for ( jj = ii-1; jj >= 0 ; jj--) {
				if (Node_Coordinate[ii][0]== Node_Coordinate[jj][0] && Node_Coordinate[ii][1]==Node_Coordinate[jj][1]) {
					//printf("同じ座標の番号ii:%d jj:%d\n",ii,jj);
					same_point[ii]=jj;
					//printf("same_point_1[%d]:%d\n",ii,same_point[ii]);
				}
			}
		}

		for (ii = 0; ii < *Total_Control_Point; ii++) {
			printf("same_point[%d]:%d\n",ii,same_point[ii]);
		}*/
//INC\の計算（節点番号をξ、ηの番号で表す為の配列）
	if (DIMENSION == 2) {
			e=0;
			for ( l = 0; l < *No_Patch; l++) {
				i=0;
			for (jj = 0; jj < No_Control_point[l][1]; jj++){
				for (ii = 0; ii < No_Control_point[l][0]; ii++) {

					INC[l][Patch_controlpoint[l][i]][0] = ii;
					INC[l][Patch_controlpoint[l][i]][1] = jj;
					//printf("INC[%d][0]=%d INC[%d][1]=%d\n",i, INC[i][0], i,INC[i][1] );
                    Adress_Controlpoint[l][ii][jj]=Patch_controlpoint[l][i];

					if (ii>= Order[l][0] && jj >= Order[l][1]) {

						for (jjloc = 0; jjloc <= Order[l][1]; jjloc++) {
							for (iiloc = 0; iiloc <= Order[l][0]; iiloc++) {
								//printf("jjloc:%d iiloc:%d\n",jjloc,iiloc);
								B = Patch_controlpoint[l][i-jjloc*No_Control_point[l][0]-iiloc];
								b = jjloc*(Order[l][0]+1)+iiloc;
								//printf("B=%d b=%d e=%d\n",B,b,e);
								Controlpoint_of_Element_nomerge[e][b] = B;
								Controlpoint_of_Element[e][b] = B;

							}
						}
						Element_patch[e]=l;
						e++;
					}
					i++;
				}
			}
		}
			/*for ( i = 0; i < *Total_Control_Point; i++){
				for ( l = 0; l < *No_Patch; l++) {
					printf("INC[%d][%d][0]=%d INC[%d][%d][1]=%d\n",l,i, INC[l][i][0], l,i,INC[l][i][1] );
				}
		}*/

		//for ( l = 0; l < *No_Patch; l++) {
		/*	for (i = 0; i < *Total_Element; i++) {
					for (b = 0; b < No_Control_point_ON_ELEMENT[Element_patch[i]]; b++) {
					printf("Controlpoint_of_Element[%d][%d]=%d\n",i,b,Controlpoint_of_Element[i][b]);
				}
			}*/
		//}

		/*for ( i = 0; i < *Total_Element; i++) {
			printf("Element_patch[%d]:%d\n",i,Element_patch[i]);
		}*/
			/*for (i = 0; i < *Total_Element; i++) {
				for (b = 0; b < No_Control_point_ON_ELEMENT[Element_patch[i]]; b++) {
					//printf("Controlpoint_of_Element_before[%d][%d]=%d\n",i,b,Controlpoint_of_Element[i][b]);
					Controlpoint_of_Element[i][b]=same_point[Controlpoint_of_Element[i][b]];
					//printf("Controlpoint_of_Element_after[%d][%d]=%d\n",i,b,Controlpoint_of_Element[i][b]);
				}
			}*/
		}

		if (DIMENSION == 3) {
			e=0;
			for (l = 0; l < *No_Patch; l++) {
				i=0;
			for (kk = 0; kk < No_Control_point[l][2] ; kk++) {
				for (jj = 0; jj < No_Control_point[l][1]; jj++){
					for (ii = 0; ii < No_Control_point[l][0]; ii++) {

						//printf("kk=%d\n",kk );
						INC[l][Patch_controlpoint[l][i]][0] = ii;
						INC[l][Patch_controlpoint[l][i]][1] = jj;
						INC[l][Patch_controlpoint[l][i]][2] = kk;
						//printf("INC[%d][0]=%d INC[%d][1]=%d\n",i, INC[i][0], i,INC[i][1] );
						if (ii>= Order[l][0] && jj >= Order[l][1] && kk >= Order[l][2]) {
							for (kkloc = 0; kkloc < Order[l][2]; kkloc++) {
								for (jjloc = 0; jjloc <= Order[l][1]; jjloc++) {
									for (iiloc = 0; iiloc <= Order[l][0]; iiloc++) {
										//printf("jjloc:%d iiloc:%d\n",jjloc,iiloc);
										B = Patch_controlpoint[l][i-jjloc*No_Control_point[l][0]-iiloc];
										b = jjloc*(Order[l][0]+1)+iiloc;
										//printf("B=%d b=%d e=%d\n",B,b,e);
										Controlpoint_of_Element_nomerge[e][b] = B;
										Controlpoint_of_Element[e][b] = B;

									}
								}
							}
							Element_patch[e]=l;
									e++;
						}
						i++;
					}
				}
			}
		}
	}
/*----------------------------------------------------------------------------------------------*/

//#include<stdio.h>

//#define DIMENSION               2
//#define MAX_N_KNOT 				1000
//#define MAX_N_ELEMENT 		    110000

//int main()
//{
    //static double Position_Knots[DIMENSION][MAX_N_KNOT];    /*ノットベクトル*/
    //static double difference[MAX_N_KNOT][DIMENSION];    /*隣り合うノットベクトルの差*/
    //static int ENC[MAX_N_ELEMENT][DIMENSION];   /*ENC[全ての要素][0,1]=x,y方向の何番目の要素か*/
    //int /*i,j,k,e,*/m,n,h,p,q,x,y;
    //int kk,rr;
    //int Order[DIMENSION];   /*次数*/
    //int No_knot[DIMENSION];   /*ノット数*/
    //int *Total_Element;  /*ゼロエレメントを含むすべての要素の数*/
    //int Total_element_all_ID[MAX_N_ELEMENT];    /*ゼロエレメントではない要素＝１、ゼロエレメント＝０*/
    //int line_No_Total_element[DIMENSION];   /*ゼロエレメントを含むすべての要素列の数*/
    //int line_No_real_element[DIMENSION];   /*ゼロエレメントではない要素列の数*/
    //int real_element_line[MAX_N_ELEMENT][DIMENSION];   /*ゼロエレメントではない要素列*/
    //int real_element[MAX_N_ELEMENT];    /*ゼロエレメントではない要素の番号*/
    //int r=0;

		for ( l = 0; l < *No_Patch; l++) {
    	for(j=0;j<DIMENSION;j++){

        //printf("Order[%d]= ",j);
        //scanf("%d",&Order[j]);

        //printf("No_knot[%d]= ",j);
        //scanf("%d",&No_knot[j]);

        line_No_Total_element[l][j]=No_knot[l][j]-2*Order[l][j]-1;

        //for(i=0;i<No_knot[j];i++){
            //printf("Position_Knots[%d][%d]= ",j,i);
            //scanf("%lf",&Position_Knots[j][i]);
        //}

        for( kkk = Order[l][j] ; kkk < No_knot[l][j]-Order[l][j]-1  ; kkk++ ){
            difference[l][kkk-Order[l][j]][j]=Position_Knots[l][j][kkk+1]-Position_Knots[l][j][kkk];
            //printf("[[[%d]]] ξ[%d]-ξ[%d]=%lf\n",kkk-Order[l][j],kkk+1,kkk,difference[kkk-Order[l][j]][j]);

            if(difference[l][kkk-Order[l][j]][j]!=0){
                line_No_real_element[l][j]++;
            }
        }
      //  printf("line_No_real_element[%d][%d]=%d\n",l,j,line_No_real_element[l][j]);
    }
	}

    /*要素に行番号、列番号をつける*/

    if(DIMENSION==2){
        //printf("Total_element_all= ");
        //scanf("%d",&Total_element_all);

        for(h=0;h<*Total_Element;h++){
            Total_element_all_ID[h]=0;
            //Total_element_all_ID[h]=h;
        }

        i=0;
				for (l = 0; l < *No_Patch; l++) {
        for(y=0;y<line_No_Total_element[l][1];y++){
            for(x=0;x<line_No_Total_element[l][0];x++){
                ENC[i][0]=x;
                ENC[i][1]=y;
                i++;

            }
        }
			}
			//int m;
        /*for(m=0;m<*Total_Element;m++){
            printf("ENC[%d][0]=%d ENC[%d][1]=%d\n",m,ENC[m][0],m,ENC[m][1]);
        }*/
    }

    /*必要な要素の行と列の番号を求める*/

    for(j=0;j<DIMENSION;j++){

				for (l = 0; l < *No_Patch; l++) {
					e=0;
					for(k=0;k<line_No_Total_element[l][j];k++){
						//printf("//%d,%d//\n",j,line_No_Total_element[j]);
						//printf("difference[%d][%d][%d]:%le\n",l,k,j,difference[l][k][j]);
						if(difference[l][k][j]!=0){
							//printf("k=%d e=%d\n",k,e);
							real_element_line[l][e][j]=k;
							e++;

						}
					}



        printf("line_No_real_element[%d][%d]:%d\n",l,j,line_No_real_element[l][j]);

       /*for(e=0;e<line_No_real_element[l][j];e++){
            printf("real_element_line[%d][%d][%d]=%d\n",l,e,j,real_element_line[l][e][j]);
        }*/
			}
    }
		/*for ( e = 0; e < *Total_Element; e++) {
			for ( j = 0; j < DIMENSION; j++) {
			printf("real_element_line[%d][%d]=%d\n",e,j,real_element_line[e][j]);
			}
		}*/
    /*必要な要素列上の要素のIDを1にする*/

    if(DIMENSION==2){
        for(n=0;n<*Total_Element;n++){
            //for(j=0;j<DIMENSION;j++){
            for(p=0;p<line_No_real_element[Element_patch[n]][0];p++){
                if(ENC[n][0]==real_element_line[Element_patch[n]][p][0]){
                    for(q=0;q<line_No_real_element[Element_patch[n]][1];q++){
                        if(ENC[n][1]==real_element_line[Element_patch[n]][q][1]){
                            //Total_element_all_ID[n]++;
                            Total_element_all_ID[n]++;
                            //break;
                        }
                    }
                }
            }
            //}

            //printf("Total_element_all_ID[%d]=%d\n",n,Total_element_all_ID[n]);

			/*IDが1の要素に番号を振る*/

            if(Total_element_all_ID[n]==1){
                real_element[r]=n;
                //printf("real_[%d]=%d\n",r,real_element[r]);
                r++;
            }
        }
				for ( l = 0; l < *No_Patch; l++) {
					real_Total_Element+=line_No_real_element[l][0]*line_No_real_element[l][1];
				}
    }
		//int rr;
 printf("real_Total_Element:%d\n",real_Total_Element);
   /*for(rr=0;rr<real_Total_Element;rr++){
        printf("real_element[%d]=%d\n",rr,real_element[rr]);
    }*/
//}
//
/* For distributed load 2D */

      for(iii = 0; iii <  *Total_Control_Point; iii++){
        Equivalent_Nodal_Force[iii][0] = 0.0;
        Equivalent_Nodal_Force[iii][1] = 0.0;
      }


    for(i=0; i < *Total_DistributeForce; i++){

      type_load = type_load_array[i];
      iPatch = iPatch_array[i];
      iCoord = iCoord_array[i];
      val_Coord = val_Coord_array[i];
      Range_Coord[0] = Range_Coord_array[i][0];
      Range_Coord[1] = Range_Coord_array[i][1];
      Coeff_Dist_Load [0] = Coeff_Dist_Load_array[i][0];
      Coeff_Dist_Load [1] = Coeff_Dist_Load_array[i][1];
      Coeff_Dist_Load [2] = Coeff_Dist_Load_array[i][2];

      Setting_Dist_Load_2D(*Total_Control_Point, iPatch, *Total_Element, iCoord, val_Coord,
			   Range_Coord, type_load, Coeff_Dist_Load);
    }
/*-------------------------------------------------------------------------------------*/

}

//拘束されている行数を省いた行列の番号の制作

int Make_Index_Dof( int Total_Control_Point, int Total_Constraint, int Constraint_Node_Dir[MAX_N_CONSTRAINT][2] ){
	int i,k=0;
	//拘束されている自由度(Degree Of free)をERRORにする
	for( i = 0; i < Total_Constraint; i++ )
		Index_Dof[ Constraint_Node_Dir[i][0]*DIMENSION + Constraint_Node_Dir[i][1] ] = ERROR;
	//ERROR以外に番号を付ける
	for( i= 0; i < Total_Control_Point*DIMENSION; i++){
		if( Index_Dof[i] != ERROR ){
			Index_Dof[i] = k;
			k++;
		}
	}
	printf("Max_Index_Dof=%d\n",k);
	return k;
}

void Make_K_Whole_Ptr_Col(int Total_Element, int Total_Control_Point, int K_Whole_Size ){
	int i,ii,j,jj,k;
	int NE;
	int N,i_index,j_index;

	//初期化
	for(i = 0; i < Total_Control_Point*DIMENSION; i ++ )
		Total_Control_Point_To_Node[i] = 0;
	for(i =0; i< K_Whole_Size+1; i++ )
		K_Whole_Ptr[ i ] = 0;

	for( N = 0; N < Total_Control_Point; N += K_DIVISION_LENGE  ){			//大きく分割するためのループ 結局N=0?
		//各節点に接する節点を取得
		for (i = 0; i < K_DIVISION_LENGE; i++)
		{
			Total_Control_Point_To_Node[i] = 0;
		}

		for( i = 0; i < Total_Element; i++ ){
			for( ii = 0; ii < No_Control_point_ON_ELEMENT[Element_patch[i]]; ii++ ){
				NE = Controlpoint_of_Element[i][ii] - N ;
				//printf("K_DIVISION_LENGE=%d,N=%d,NE=%d\n",K_DIVISION_LENGE,N,NE);    //K_DIVISION_LENGE=0,N=0,NE=コネクティビティ的な
				if( 0 <= NE && NE < K_DIVISION_LENGE ){
					for( j = 0; j <No_Control_point_ON_ELEMENT[Element_patch[i]]; j++ ){
						//printf("j=%d\n",j);
						//数字がない時
						if( Total_Control_Point_To_Node[ NE ] == 0 ){
							//節点番号を取得
							Node_To_Node[ NE ][ 0 ] = Controlpoint_of_Element[i][j];
							Total_Control_Point_To_Node[ NE ] ++ ;
							//printf("i=%d,ii=%d,j=%d,Total_Control_Point_To_Node[ NE ] == 0のif文やった\n",i,ii,j);
							//printf("Node_To_Node[%d][0]=%d\n",NE,Node_To_Node[NE][0]);
							//printf("①Total_Control_Point_To_Node[%d]=%d\n",NE,Total_Control_Point_To_Node[NE]);
						}
						//printf("②Total_Control_Point_To_Node[%d]=%d\n",NE,Total_Control_Point_To_Node[NE]);
						//同じものがあったら
						for( k = 0; k < Total_Control_Point_To_Node[ NE ]; k++ )
							//printf("Node_To_Node[%d][%d]=%d\n",NE,k,Node_To_Node[NE][k]);絶対コメントアウト　出力結果変わっちゃう
							if( Node_To_Node[ NE ][k] == Controlpoint_of_Element[i][j] )break;
						if( k == Total_Control_Point_To_Node[ NE ] ){
							Node_To_Node[ NE ][k] = Controlpoint_of_Element[i][j];
							Total_Control_Point_To_Node[ NE ] ++ ;
							//printf("i=%d,ii=%d,j=%d,k=%d,k=Total_Control_Point_To_Node[ NE ]のif文やった\n",i,ii,j,k);
							//printf("③Total_Control_Point_To_Node[%d]=%d\n",NE,Total_Control_Point_To_Node[NE]);
						}
					}
				}
			}
			//printf("\n");
		}
		//順番に並び替える
		for( i = 0; i < K_DIVISION_LENGE; i++ ){
			if( N+i < Total_Control_Point ){
				//printf("Node[%d] T=%d; \n",N+i, Total_Control_Point_To_Node[ i ]);
				for( j = 0; j < Total_Control_Point_To_Node[ i ]; j++){
					int Min = Node_To_Node[i][j], No=j;
					for( k = j; k < Total_Control_Point_To_Node[ i ]; k++ ){
						if( Min >  Node_To_Node[i][k] ){
							Min = Node_To_Node[i][k];	No = k;
						}
					}
					for( k = No; k > j; k-- )
						Node_To_Node[i][k] = Node_To_Node[i][k-1];
					Node_To_Node[i][j] = Min;
	//				printf("%d ",Node_To_Node[i][j]);
				}
	//			printf("\n");
			}
		}


		//節点からcol ptrを求める
		ii = 0; k= 0;
		for( i = 0; i < K_DIVISION_LENGE; i++ ){
			for( ii = 0; ii < DIMENSION; ii++ ){
				if( N+i < Total_Control_Point ){
					i_index = Index_Dof[(N+i)*DIMENSION+ii];
					k = 0;
					if( i_index >= 0 ){
						K_Whole_Ptr[ i_index + 1 ] = K_Whole_Ptr[ i_index ];
						for( j = 0; j < Total_Control_Point_To_Node[i]; j++ ){
							for(jj = 0; jj < DIMENSION; jj++ ){
								j_index = Index_Dof[ Node_To_Node[i][j]*DIMENSION+jj ];
								if( j_index >= 0 && j_index >= i_index ){
									K_Whole_Ptr[ i_index + 1 ] ++;
									//col_N[N/K_DIVISION_LENGE][k] = j_index;
									K_Whole_Col[ K_Whole_Ptr[ i_index ] + k ] = j_index;
									k++;
									//printf("ptr[%d]=%d,col[%d]=%d\n",i_index+1,K_Whole_Ptr[i_index+1],K_Whole_Ptr[i_index]+k,K_Whole_Col[K_Whole_Ptr[i_index]+k]);

								}
							}
						}
					}
				}
			}
		}
		//col_N[N/K_DIVISION_LENGE][ k ] = -1;
	}
	/*
	for( i = 0; i < K_Whole_Size+1; i++ )//printf("K_Whole_Ptr[%d]= %d\n",i,K_Whole_Ptr[i]);
	//col合成
	k = 0;
	for( N = 0; N < Total_Control_Point ; N +=K_DIVISION_LENGE ){
		for(i = 0; col_N[ N/K_DIVISION_LENGE ][i] != -1; i++ ){
			K_Whole_Col[k] = col_N[ N/K_DIVISION_LENGE ][i];
			k++;
		}
	}*/
}
//valを求める
void Make_K_Whole_Val( double E, double nu, int Total_Element, int K_Whole_Size ,int DM,int Total_Control_Point/*,int real_element[MAX_N_ELEMENT]*/){
	int i, j1, j2, k1, k2, l;
		int a,b,re;

		for( i = 0; i < MAX_NON_ZERO; i++ ) K_Whole_Val[i] = 0.0;

		/*for(rr=0;rr<line_No_real_element[0]*line_No_real_element[1];rr++){
        	printf("real_element[%d]=%d\n",rr,real_element[rr]);
    	}
		printf("real_Total_element=%d\n",real_Total_Element);*/
		//re=0;
		for(re = 0; re < real_Total_Element; re++ ){

			//if(i==real_element[re]){
				//printf("re=%d\n",re);
			i=real_element[re];
			KIEL_SIZE=No_Control_point_ON_ELEMENT[Element_patch[i]]*DIMENSION;
			double X[No_Control_point_ON_ELEMENT[Element_patch[i]]][DIMENSION], K_EL[KIEL_SIZE][KIEL_SIZE];
			//printf("El_No=%d\n",i );
			//各要素のKelを求める
			for( j1 = 0; j1 < No_Control_point_ON_ELEMENT[Element_patch[i]]; j1++ ){
				for( j2 = 0; j2 < DIMENSION; j2++ ){
					X[j1][j2] = Node_Coordinate[ Controlpoint_of_Element[i][j1] ][ j2 ];
				}
			}

			Make_K_EL( i, X, K_EL, E, nu ,DM , Total_Element,Total_Control_Point);

			//Valを求める
			for( j1 = 0; j1 < No_Control_point_ON_ELEMENT[Element_patch[i]]; j1++ ){for( j2 = 0; j2 < DIMENSION; j2++ ){
				a = Index_Dof[ Controlpoint_of_Element[i][j1]*DIMENSION+j2 ];
				if( a >= 0 ){
					for( k1 = 0; k1 < No_Control_point_ON_ELEMENT[Element_patch[i]]; k1++ ){for( k2 = 0; k2 < DIMENSION; k2++ ){
						b = Index_Dof[ Controlpoint_of_Element[i][k1]*DIMENSION+k2 ];
						if( b >= 0 && b >= a ){
							for( l = K_Whole_Ptr[a]; l <  K_Whole_Ptr[a+1]; l++){
								if( K_Whole_Col[l] == b ){
									K_Whole_Val[l] += K_EL[j1*DIMENSION+j2][k1*DIMENSION+k2];
									//printf("l;%d\n",l);
									break;
								}
							}
						}
					}}
				}
			}}

			//re++;
			//}

		}
}

///////////////////////////////////////////////////////////////////////////
/////////////////////連立1次方程式の解法
/////////////////////////////////////////////////////////////////////
//分布荷重の等価節点力を足す
void Add_Equivalent_Nodal_Forec_to_F_Vec(int Total_Control_Point){
  int i,j,index;
  for(j=0; j < DIMENSION; j++)
	for( i = 0; i < Total_Control_Point; i++ ){
		index = Index_Dof[i * DIMENSION + j];
		if( index >= 0 ){
			rhs_vec[index] += Equivalent_Nodal_Force[i][j];
			 //printf("i = %d index = %d rhs_vec[index] = %f\n",i, index, rhs_vec[index]);
			}
			//printf("i = %d  j = %d  Equivalent_Nodal_Force[i][j] = %f\n",i, j, Equivalent_Nodal_Force[i][j]);
	}
}
//荷重の行列を作る
void Make_F_Vec( int Total_Load,int Load_Node_Dir[MAX_N_LOAD][2],double Value_of_Load[MAX_N_LOAD],int K_Whole_Size ){
	int i,index;
	for( i = 0; i < K_Whole_Size; i++ )
		rhs_vec[i] = 0.0;
	for( i = 0; i < Total_Load; i++ ){
		index = Index_Dof[ Load_Node_Dir[i][0] * DIMENSION + Load_Node_Dir[i][1] ];
		if( index >= 0 )
			rhs_vec[index] += Value_of_Load[i];
	}
}
//強制変位対策
void Make_F_Vec_disp_const(int Total_Constraint,int Constraint_Node_Dir[MAX_N_CONSTRAINT][2],double Value_of_Constraint[MAX_N_CONSTRAINT],int Total_Element,double  E,double nu, int DM, int Total_Control_Point){
		 int ie , idir, inode, jdir, jnode, kk_const;
		 int ii,iii,b,bb,jj,j1,j2,ii_local,jj_local;

		 double K_EL[KIEL_SIZE][KIEL_SIZE];

			for(ie=0; ie < Total_Element; ie++){
				double X[No_Control_point_ON_ELEMENT[Element_patch[ie]]][DIMENSION];
				iii=0;
				for(idir=0; idir < DIMENSION ; idir++){
					for(inode=0; inode < No_Control_point_ON_ELEMENT[Element_patch[ie]]; inode++){
						b = Index_Dof[ Controlpoint_of_Element[ie][inode]*DIMENSION+idir ];
						if(b < 0) iii++;
					}
				}
				if(iii > 0){
					for( j1 = 0; j1 < No_Control_point_ON_ELEMENT[Element_patch[ie]]; j1++ ){
						for( j2 = 0; j2 < DIMENSION; j2++ ){
							X[j1][j2] = Node_Coordinate[ Controlpoint_of_Element[ie][j1]][j2];
						}//end for j2
					}//end for j1
					Make_K_EL( ie, X, K_EL, E, nu , DM, Total_Element,Total_Control_Point);
					for(idir=0; idir < DIMENSION ; idir++){
						for(inode=0; inode < No_Control_point_ON_ELEMENT[Element_patch[ie]]; inode++){
							ii = Controlpoint_of_Element[ie][inode]*DIMENSION+idir;
							b = Index_Dof[ ii ];
								if(b >= 0 ){
									ii_local = inode * DIMENSION +idir;
									for(jdir=0; jdir < DIMENSION ; jdir++){
										for(jnode=0; jnode < No_Control_point_ON_ELEMENT[Element_patch[ie]]; jnode++){
											jj = Controlpoint_of_Element[ie][jnode]*DIMENSION+jdir;
											bb = Index_Dof[ jj ];
											if(bb < 0){
												jj_local = jnode * DIMENSION +jdir;//printf("%d,%d\n",ie,jnode);
												for(kk_const = 0; kk_const < Total_Constraint ;  kk_const++){
													if(Controlpoint_of_Element[ie][jnode] == Constraint_Node_Dir[kk_const][0] && jdir == Constraint_Node_Dir[kk_const][1]){
														rhs_vec[b] -= K_EL[ii_local][jj_local]* Value_of_Constraint[kk_const];//if(kk_const >= 28){printf("%d , %d ,%16.15e\n",ii_local, jj_local ,  K_EL[ii_local][jj_local]);}
													}//end if Controlpoint_of_Element[ie][jnode]
												}//end for kk_const
											}//end if bb
										}//end for jnode
									}//end for jdir
								}//end if b>=0
							}//end for inode
						}//end for idir
					}//end if iii>0
				}//end for ie
}//end

void mat_vec_crs(double vec_result[], double vec[], const int ndof){
	int i,j,icount = 0;
	/* zero clear */

	for(i=0; i < ndof; i++) vec_result[i] = 0;
	for(i=0; i < ndof; i++){
		for(j=K_Whole_Ptr[i]; j < K_Whole_Ptr[i+1]; j++){
			vec_result[i] += K_Whole_Val[icount] * vec[K_Whole_Col[j]];
			if(i != K_Whole_Col[j]) vec_result[K_Whole_Col[j]] += K_Whole_Val[icount] * vec[i];
			icount ++;

		}
	}
}
double inner_product(int ndof, double vec1[], double vec2[]){
	double rrr=0.0;
	int i;
	for(i=0; i<ndof; i++){
         rrr += vec1[i]*vec2[i];
         //printf("vec1[%d]=%f vec2[%d]=%f\n",i,vec1[i],i,vec2[i]); /*-nan 10/23*/
   }
	return(rrr);

}
int check_conv_CG(int ndof, double alphak, double pp[], double eps, int itr){
	double rrr1=0.0, rrr2=0.0, rrr3;
	int i, istop=0;
	/* Checking the convergence of the CG solver */
	/* istop =0; Not converged, istop = 1; converged */
	printf("ndof=%d alphak= %15e\t",ndof,alphak );
	for(i=0; i < ndof; i++){
		rrr1 += pp[i]*pp[i];
		rrr2 += sol_vec[i]*sol_vec[i];
        //printf("pp[%d]=%f sol_vec[%d]=%f\n",i,pp[i],i,sol_vec[i]); /*-nan 10/23*/
	}
	rrr3 = fabs(alphak) * sqrt(rrr1/rrr2);
	printf("Iteration# = %d  residual = %15e (%15e)\n",itr, rrr3, eps);
	if(rrr3 < eps) istop = 1;
	/* Temporaty Oct. 10, 2017 by H.Okada */
	//if(itr < 100) istop=0;
	//printf("Iteration# = %d  residual = %15e (%15e)\n",itr, rrr3, eps);
	return(istop);
}
void Diag_Scaling_CG_pre(int ndof, int flag_operation){
	int i,j;
	int icount = 0;
	/* flag_opertion = 0: Preprocess to the CG solver
			A <-- Dt A D  and b <-- Dt b */
	/* flag_operation = 1: Post process to the CG solver
			b <-- Dt b  */
	if(flag_operation == 0){
		diag_scaling[0] = 1.0/sqrt(K_Whole_Val[0]);
		/* diag_scaling[0] = 1.0; */
		for(i=1; i<ndof; i++){
			//printf("%d %le\n",K_Whole_Ptr[i], K_Whole_Val[K_Whole_Ptr[i]]);
			diag_scaling[i] = 1.0/sqrt(K_Whole_Val[K_Whole_Ptr[i]]);
			/* diag_scaling[i] = 1.0; */
		}
		for(i=0; i < ndof; i++){
			for(j=K_Whole_Ptr[i]; j < K_Whole_Ptr[i+1]; j++){
				//printf("Check scling icount=%d i=%d K_Whole_Col[%d] = %d\n",icount,i,j,K_Whole_Col[j]);
				K_Whole_Val[icount] = K_Whole_Val[icount]* diag_scaling[i] * diag_scaling[K_Whole_Col[j]];
				//printf("K_Whole_Val = %f\n",K_Whole_Val[icount]);
				icount ++;
			}
			//printf("rhs_vec_before[%d]:%le diag_scaling[%d]:%le\n", i,rhs_vec[i],i,diag_scaling[i]);
			rhs_vec[i] = rhs_vec[i] * diag_scaling[i];
			//printf("rhs_vec[%d]:%le\n",i,rhs_vec[i]);
		}
	}
	if(flag_operation == 1)
		for(i=0; i < ndof; i++){
			//printf("solvec[%d] = %f\n",i, sol_vec[i]);
			sol_vec[i] = sol_vec[i] *  diag_scaling[i];
		}

		printf("\nqq\n");
}
void CG_Solver(int ndof, int max_itr, double eps, int flag_ini_val){
	static double gg[MAX_K_WHOLE_SIZE], dd[MAX_K_WHOLE_SIZE], pp[MAX_K_WHOLE_SIZE];
	static double qqq, ppp, rrr;
	static double alphak, betak;
	int i;
	int itr;
	int ii, istop;

	/* Program to solve linear equations by using the CG method */
	if(flag_ini_val ==0)
    for(i=0; i < ndof; i++) sol_vec[i] = 0.0;
	/* Initializing the solution vector if it were not given */
	mat_vec_crs(dd, sol_vec, ndof);
	for(i=0; i < ndof; i++){
		gg[i] = rhs_vec[i] - dd[i];
		//printf("rhs_vec[%d]=%f dd[%d]=%f\n",i,rhs_vec[i],i,dd[i]);
		pp[i] = gg[i];
	}

			printf("\nrr");

	for(itr=0; itr < max_itr; itr++){
		ppp = inner_product(ndof, gg, gg);
		mat_vec_crs(dd, pp, ndof);
		rrr = inner_product(ndof, dd, pp);
		alphak = ppp / rrr;
        //printf("ppp=%f rrr=%f\n", ppp, rrr); /*ppp,rrrも-nan,10/22*/
        //printf("i=%d",i);
		for(ii = 0; ii < ndof; ii++){
			sol_vec[ii] += alphak * pp[ii];
			gg[ii] -= alphak * dd[ii];
		}
		qqq = inner_product(ndof, gg, dd);
		betak = qqq / rrr;
		for(ii = 0; ii < ndof; ii++) pp[ii] = gg[ii] - betak * pp[ii];
		istop = check_conv_CG(ndof, alphak, pp, eps, itr);
		if(istop == 1) break;
	}

			printf("\nss");
}
////////////////////////////////////////////////////////////////////////
/////////////////要素剛性マトリックス
////////////////////////////////////////////////////////////////////////
//IGAの基底関数
void ShapeFunction1D(double Position_Data_param[DIMENSION], int j,int e){

	int ii;
	int p;

	//printf("shapefuc_Position_Data_param[%d]:%le\n", j,Position_Data_param[j]);

	for(ii=0; ii<No_knot[Element_patch[e]][j];ii++){
		if (Position_Knots[Element_patch[e]][j][ii]==Position_Knots[Element_patch[e]][j][ii+1]) {
			Shape[j][ii][0] = 0.0;
		}
		else if (Position_Knots[Element_patch[e]][j][ii] != Position_Knots[Element_patch[e]][j][ii+1] && Position_Knots[Element_patch[e]][j][ii]<=Position_Data_param[j] && Position_Data_param[j]< Position_Knots[Element_patch[e]][j][ii+1]){
			Shape[j][ii][0] = 1.0;
		}
		else if (Position_Knots[Element_patch[e]][j][ii] != Position_Knots[Element_patch[e]][j][ii+1] && Position_Knots[Element_patch[e]][j][ii+1]== Position_Knots[Element_patch[e]][j][(No_knot[Element_patch[e]][j]-1)] &&Position_Knots[Element_patch[e]][j][ii]<=Position_Data_param[j] && Position_Data_param[j]<= Position_Knots[Element_patch[e]][j][ii+1]) {
				Shape[j][ii][0] = 1.0;
		}
		else Shape[j][ii][0]=0.0;
		//printf("Shape[%d][%d][0]=%le   ",j,ii,Shape[j][ii][0]);
	}

			for(ii=0; ii<No_knot[Element_patch[e]][j];ii++){
				for (p = 1; p <=Order[Element_patch[e]][j]; p++) {
					Shape[j][ii][p] = 0;

			}
		}
		double left_term, right_term;
		for(p=1; p<=Order[Element_patch[e]][j]; p++){
		 for(ii=0; ii<No_knot[Element_patch[e]][j]; ii++){
			 left_term = 0.0;
			 right_term = 0.0;

			 if((Position_Data_param[j]-Position_Knots[Element_patch[e]][j][ii])*Shape[j][ii][p-1] == 0 && Position_Knots[Element_patch[e]][j][ii+p]-Position_Knots[Element_patch[e]][j][ii] == 0) left_term = 0.0;
			 else left_term = (Position_Data_param[j]-Position_Knots[Element_patch[e]][j][ii])/(Position_Knots[Element_patch[e]][j][ii+p]-Position_Knots[Element_patch[e]][j][ii])*Shape[j][ii][p-1];

			 if((Position_Knots[Element_patch[e]][j][ii+p+1]-Position_Data_param[j])*Shape[j][ii+1][p-1] == 0 && Position_Knots[Element_patch[e]][j][ii+p+1]-Position_Knots[Element_patch[e]][j][ii+1] == 0) right_term = 0.0;
			 else right_term = (Position_Knots[Element_patch[e]][j][ii+p+1]-Position_Data_param[j])/(Position_Knots[Element_patch[e]][j][ii+p+1]-Position_Knots[Element_patch[e]][j][ii+1])*Shape[j][ii+1][p-1];

						Shape[j][ii][p] = left_term + right_term;
						//printf("Shape[%d][%d][%d]=%le   ",j,ii,p,Shape[j][ii][p]);
					}
		}
		//printf("order[%d]:%d\n",j,Order[Element_patch[e]][j] );
		double dleft_term,dright_term;
			for (ii = 0; ii < No_Control_point[Element_patch[e]][j]+1; ii++) {
				//printf("No_Control_point[%d]=%d\n",j,No_Control_point[j] );
		dleft_term=0.0;
		dright_term=0.0;

		if (Order[Element_patch[e]][j]*Shape[j][ii][Order[Element_patch[e]][j]-1] == 0 && Position_Knots[Element_patch[e]][j][ii+Order[Element_patch[e]][j]]-Position_Knots[Element_patch[e]][j][ii] == 0) dleft_term = 0.0 ;
		else dleft_term =  Order[Element_patch[e]][j]/(Position_Knots[Element_patch[e]][j][ii+Order[Element_patch[e]][j]]-Position_Knots[Element_patch[e]][j][ii])*Shape[j][ii][Order[Element_patch[e]][j]-1];
		/*printf("test_Shape_left[%d][%d][%d]=%le\n", j,ii,Order[Element_patch[e]][j]-1,Shape[j][ii][Order[Element_patch[e]][j]-1]);
		printf("Position_Knots[Element_patch[e]][%d][%d]:%le\n", j,ii+Order[Element_patch[e]][j],Position_Knots[Element_patch[e]][j][ii+Order[Element_patch[e]][j]]);
		printf("Position_Knots[Element_patch[e]][%d][%d]:%le\n", j,ii,Position_Knots[Element_patch[e]][j][ii]);
		printf("dleft_term=%f\n",dleft_term );*/

		if (Order[Element_patch[e]][j]*Shape[j][ii+1][Order[Element_patch[e]][j]-1] == 0 && Position_Knots[Element_patch[e]][j][ii+Order[Element_patch[e]][j]+1]-Position_Knots[Element_patch[e]][j][ii+1] == 0) dright_term = 0.0 ;
		else dright_term = Order[Element_patch[e]][j]/(Position_Knots[Element_patch[e]][j][ii+Order[Element_patch[e]][j]+1]-Position_Knots[Element_patch[e]][j][ii+1])*Shape[j][ii+1][Order[Element_patch[e]][j]-1];
		/*printf("test_Shape_right[%d][%d][%d]=%le\n", j,ii+1,Order[Element_patch[e]][j]-1,Shape[j][ii+1][Order[Element_patch[e]][j]-1]);
		printf("Position_Knots[%d][%d]:%le\n", j,ii+Order[Element_patch[e]][j]+1,Position_Knots[j][ii+Order[Element_patch[e]][j]+1]);
		printf("Position_Knots[%d][%d]:%le\n", j,ii+1,Position_Knots[j][ii+1]);
		printf("dright_term=%f\n",dright_term );*/

		dShape[j][ii]= dleft_term - dright_term;

		//printf("PP=%d\n",PP );

		//printf("dShape[%d][%d]= %f\n",j,ii,dShape[j][ii]);
	}



}

void ShapeFunc_from_paren(double Local_coord[DIMENSION], int j, int e){

	int i=0;
  //printf("Local_coord[%d]:%le\n",j,Local_coord[j]);
	i=INC[Element_patch[e]][Controlpoint_of_Element_nomerge[e][0]][j];
	//printf("El_No:%d\n",e );
	//printf("i:%d\n",i);
	//printf("Position_Knots[%d][%d]:%le Position_Knots[%d][%d]:%le\n", j,i+1,Position_Knots[j][i+1],j,i,Position_Knots[j][i]);
	Position_Data_param[j] = ((Position_Knots[Element_patch[e]][j][i+1]-Position_Knots[Element_patch[e]][j][i])* Local_coord[j]+(Position_Knots[Element_patch[e]][j][i+1]+Position_Knots[Element_patch[e]][j][i]))/2;
	//printf("Position_Data_param[%d]:%le\n",j,Position_Data_param[j]);
}

double dShapeFunc_from_paren(int j,int e) {
	int i;
	double dPosition_Data_param;

	i=INC[Element_patch[e]][Controlpoint_of_Element_nomerge[e][0]][j];
	//printf("El_No:%d\n",e );
	//printf("i:%d\n",i);
	//printf("Position_Knots[Element_patch[e]][%d][%d]:%le Position_Knots[%d][%d]:%le\n", j,i+1,Position_Knots[j][i+1],j,i,Position_Knots[j][i]);
	dPosition_Data_param = (Position_Knots[Element_patch[e]][j][i+1]-Position_Knots[Element_patch[e]][j][i])/2;
	//printf("dPosition_Data_param:%le\t",dPosition_Data_param);
	return dPosition_Data_param;
}


double Shape_func (int I_No,int Total_Control_Point, double Local_coord[DIMENSION],int El_No) {
//printf("El_No:%d I_No:%d\n",El_No,I_No);
	int i,j;
	double R;
	double weight_func;
		weight_func=0.0;
		//shape_func[]={0.0};

		for (i = 0; i < Total_Control_Point; i++) {
			shape_func[i]=1.0;
		}



		for (j = 0; j < DIMENSION; j++) {
			ShapeFunc_from_paren(Local_coord,j,El_No);
			ShapeFunction1D(Position_Data_param,j,El_No);
			for ( i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++) {
				shape_func[Controlpoint_of_Element_nomerge[El_No][i]] *= Shape[j][INC[Element_patch[El_No]][Controlpoint_of_Element_nomerge[El_No][i]][j]][Order[Element_patch[El_No]][j]];    /*基底関数*/
				//printf("%d",shape_func[0][0]);
			}
		}

		for ( i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++) {
			weight_func += shape_func[Controlpoint_of_Element_nomerge[El_No][i]]*Node_Coordinate[Controlpoint_of_Element[El_No][i]][DIMENSION];
			/*if (El_No==2 && I_No==0) {
			printf("Controlpoint_of_Element_nomerge[%d][%d]:%d shape_func[%d]:%le\n", El_No,i,Controlpoint_of_Element_nomerge[El_No][i],Controlpoint_of_Element_nomerge[El_No][i],shape_func[Controlpoint_of_Element_nomerge[El_No][i]]);
			printf("Controlpoint_of_Element[%d][%d]:%d Node_Coordinate[%d][%d]:%lf\n",El_No,i, Controlpoint_of_Element[El_No][i],Controlpoint_of_Element[El_No][i],DIMENSION,Node_Coordinate[Controlpoint_of_Element[El_No][i]][DIMENSION]);
		}*/
		}
	//printf("weight_func:%le\n",weight_func );
	if(I_No < No_Control_point_ON_ELEMENT[Element_patch[El_No]]) R = shape_func[Controlpoint_of_Element_nomerge[El_No][I_No]]*Node_Coordinate[Controlpoint_of_Element[El_No][I_No]][DIMENSION]/weight_func;

	else R = ERROR;
	//printf("shape_func R:%le\n",R );
	return R;
 }

void NURBS_deriv(double Local_coord[DIMENSION],int El_No, int Total_Control_Point) {
	double weight_func;
	//double shape_func[100][50];


	double dWeight_func1;
	double dWeight_func2;

	int i,j;
	//int ii;


		//for(ii = 0; ii < NN+1; ii++){
			//printf("NdShape3[%d]= %f\n",ii,dShape3[ii]);
		//}

		//for (ii = 0; ii < NN+1; ii++)printf("NdShape1[%d]= %f\n",ii,dShape1[ii]);
		//for (jj = 0; jj < MM+1; jj++)printf("NdShape2[%d]= %f\n",jj,dShape2[jj]);
		//printf("\n");

		for (i = 0; i < Total_Control_Point; i++) {
			shape_func[i]=1.0;
		}

		weight_func=0.0;

		dWeight_func1=0.0;
		dWeight_func2=0.0;



		for (j = 0; j < DIMENSION; j++) {
			ShapeFunc_from_paren(Local_coord,j,El_No);
			ShapeFunction1D(Position_Data_param,j,El_No);
			for ( i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++) {
				shape_func[Controlpoint_of_Element_nomerge[El_No][i]]*= Shape[j][INC[Element_patch[El_No]][Controlpoint_of_Element_nomerge[El_No][i]][j]][Order[Element_patch[El_No]][j]];
				//printf("Shape[%d][%d][%d]:%le\n",j,INC[Controlpoint_of_Element[El_No][i]][j],Order[Element_patch[e]][j],Shape[j][INC[Controlpoint_of_Element[El_No][i]][j]][Order[Element_patch[e]][j]]);
			}
		}

		for ( i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++) {
			weight_func += shape_func[Controlpoint_of_Element_nomerge[El_No][i]]*Node_Coordinate[Controlpoint_of_Element_nomerge[El_No][i]][DIMENSION];
			//printf("Node_Coordinate[%d][%d]:%le\n", Controlpoint_of_Element[El_No][i],DIMENSION,Node_Coordinate[Controlpoint_of_Element[El_No][i]][DIMENSION]);
			//printf("shape_func[%d]:%le\n",Controlpoint_of_Element[El_No][i],shape_func[Controlpoint_of_Element[El_No][i]]);
			//printf("weight_func:%le\n", weight_func);
		}
	//printf("weight_func:%le\n", weight_func);
		//for(jj=0;jj<NN;jj++) for(kk=0;kk<MM;kk++)shape_func[jj][kk] = Shape1[jj][PP]*Shape2[kk][QQ]*Node_Coordinate[MM*jj+kk][2]/weight_func;
		for ( i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++){
				dWeight_func1 += dShape[0][INC[Element_patch[El_No]][Controlpoint_of_Element_nomerge[El_No][i]][0]]*Shape[1][INC[Element_patch[El_No]][Controlpoint_of_Element_nomerge[El_No][i]][1]][Order[Element_patch[El_No]][1]]*Node_Coordinate[Controlpoint_of_Element[El_No][i]][DIMENSION];
				dWeight_func2 += Shape[0][INC[Element_patch[El_No]][Controlpoint_of_Element_nomerge[El_No][i]][0]][Order[Element_patch[El_No]][0]]*dShape[1][INC[Element_patch[El_No]][Controlpoint_of_Element_nomerge[El_No][i]][1]]*Node_Coordinate[Controlpoint_of_Element[El_No][i]][DIMENSION];
		}
		//printf("dWeight_func1:%le dWeight_func2:%le\n",dWeight_func1,dWeight_func2);
		for ( i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++){
				dShape_func1[Controlpoint_of_Element_nomerge[El_No][i]]=Node_Coordinate[Controlpoint_of_Element[El_No][i]][DIMENSION]*(weight_func*dShape[0][INC[Element_patch[El_No]][Controlpoint_of_Element_nomerge[El_No][i]][0]]*Shape[1][INC[Element_patch[El_No]][Controlpoint_of_Element_nomerge[El_No][i]][1]][Order[Element_patch[El_No]][1]]-dWeight_func1*shape_func[Controlpoint_of_Element_nomerge[El_No][i]])/(weight_func*weight_func);
				dShape_func2[Controlpoint_of_Element_nomerge[El_No][i]]=Node_Coordinate[Controlpoint_of_Element[El_No][i]][DIMENSION]*(weight_func*Shape[0][INC[Element_patch[El_No]][Controlpoint_of_Element_nomerge[El_No][i]][0]][Order[Element_patch[El_No]][0]]*dShape[1][INC[Element_patch[El_No]][Controlpoint_of_Element_nomerge[El_No][i]][1]]-dWeight_func2*shape_func[Controlpoint_of_Element_nomerge[El_No][i]])/(weight_func*weight_func);
				//printf("NURBS_deriv;Controlpoint_of_Element[%d][%d]:%d\n",El_No,i,Controlpoint_of_Element[El_No][i]);
				//printf("dShape_func1[%d]:%le\n",Controlpoint_of_Element[El_No][i],dShape_func1[Controlpoint_of_Element[El_No][i]]);
				//printf("dShape_func2[%d]:%le\n",Controlpoint_of_Element[El_No][i],dShape_func2[Controlpoint_of_Element[El_No][i]]);
		}
	}


double dShape_func(int I_No, int xez, double Local_coord[DIMENSION], int El_No,int Total_Control_Point){
	double dR;


 	//printf("El_No=%d\n",El_No );

	NURBS_deriv(Local_coord,El_No,Total_Control_Point);

	if(xez!=0 && xez!=1) dR = ERROR;

		else if(I_No < No_Control_point_ON_ELEMENT[Element_patch[El_No]])
		{
			if( xez == 0 )	  {  dR = dShape_func1[Controlpoint_of_Element_nomerge[El_No][I_No]]*dShapeFunc_from_paren(xez,El_No);
		//printf("dShape_func1[%d]:%le\n",Controlpoint_of_Element[El_No][I_No],dShape_func1[Controlpoint_of_Element[El_No][I_No]]);
		}
			else if( xez == 1 ) {dR = dShape_func2[Controlpoint_of_Element_nomerge[El_No][I_No]]*dShapeFunc_from_paren(xez,El_No);
			//printf("dShape_func2[%d]:%le\n",Controlpoint_of_Element[El_No][I_No],dShape_func2[Controlpoint_of_Element[El_No][I_No]]);
		}
		}


	else dR = ERROR;
	//printf("dR:%le\n",dR);

	//printf("I_No=%d xez=%d dR=%le\t", I_No, xez, dR );
	//printf("Controlpoint_of_Element[%d][%d]:%d\t",El_No,I_No,Controlpoint_of_Element[El_No][I_No]);
	//printf("dShape_func1:%le\t",dShape_func1[Controlpoint_of_Element[El_No][I_No]]);
	//printf("dShape_func2:%le\n",dShape_func2[Controlpoint_of_Element[El_No][I_No]]);

/*for (i = 0; i < DIMENSION; i++) {
	printf("Local_coord[%d]=%lf\n",i,Local_coord[i] );
}*/
	//printf("dR:%le\t", dR);
	return dR;
}
/*
//形状関数
double N_Quad_4(int I_No, double Local_coord[DIMENSION] )
{
	double N;
	if(I_No==0) N= (1.0+Local_coord[0])*(1.0-Local_coord[1])/4.0;
	else if(I_No==1) N = (1.0+Local_coord[0])*(1.0+Local_coord[1])/4.0;
	else if(I_No==2) N = (1.0-Local_coord[0])*(1.0+Local_coord[1])/4.0;
	else if(I_No==3) N = (1.0-Local_coord[0])*(1.0-Local_coord[1])/4.0;
	else N = ERROR;
	return N;
}

//形状関数の偏微分（I_No:節点番号 xez:偏微分の分母部分0ξ1η2ζ）
double dN_Quad_4(int I_No, double Local_coord[DIMENSION], int xez)
{
	double dN;
	if(xez!=0 && xez!=1) dN = ERROR;

	else if(I_No==0)
		{
		if( xez == 0 )	    dN = (1.0-Local_coord[1])/4.0;
		else if( xez == 1 ) dN = (1.0+Local_coord[0])*(-1)/4.0;
		}

	else if(I_No==1)
		{
		if( xez == 0 )      dN = (1.0+Local_coord[1])/4.0;
		else if( xez == 1 ) dN = (1.0+Local_coord[0])/4.0;
		}

	else if(I_No==2)
		{
		if( xez == 0 )      dN = (1.0+Local_coord[1])*(-1)/4.0;
		else if( xez == 1 ) dN = (1.0-Local_coord[0])/4.0;
		}

	else if(I_No==3)
		{
		if( xez == 0 )      dN = (1.0-Local_coord[1])*(-1)/4.0;
		else if( xez == 1 ) dN = (1.0-Local_coord[0])*(-1)/4.0;
		}

	else dN = ERROR;

	return dN;
}
*/

//逆行列を元の行列に代入
double InverseMatrix_2D(double M[2][2] ){
	int i, j;
	double a[2][2];
	double det = M[0][0]*M[1][1] -M[0][1]*M[1][0];

	if(det==0) return ERROR;

	for( i= 0; i< 2; i++ )
	{
		for( j= 0; j< 2; j++ )	a[i][j] = M[i][j];
	}
	M[0][0] = a[1][1]/det;	    M[0][1] = a[0][1]*(-1)/det;
	M[1][0] = a[1][0]*(-1)/det; M[1][1] = a[0][0]/det;
	//printf("det;%le\n", det);
	return det;
}

int Jacobian( int El_No, double a[DIMENSION][DIMENSION], double Local_coord[DIMENSION], double X[No_Control_point_ON_ELEMENT[Element_patch[El_No]]][DIMENSION],
	 	 int Total_Control_Point){
	int i,j,k;
	//printf("El_No_jacobi:%d\n",El_No);
	for( i= 0; i < DIMENSION; i++ ){
		for(j = 0; j < DIMENSION; j++ ){
			a[i][j] = 0.0;
			for( k = 0; k < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; k++ ){
				//printf("Local_coord[%d]:%le\n",j,Local_coord[j] );
				a[i][j] += dShape_func(k,j,Local_coord,El_No,Total_Control_Point) * X[k][i];
				//printf(" X[%d][%d]=%lf\t",k,i, X[k][i] );
				//printf("k=%d a[%d][%d]:%le\n",k,i,j,a[i][j]);
		}
		//printf("<<<最終a[%d][%d]:%le>>>\n",i,j,a[i][j]);
		}
	}
	/*for (i = 0; i < DIMENSION; i++) {
		for (j = 0; j < DIMENSION; j++) {
			printf("a[%d][%d]:%le\n",i,j,a[i][j]);
		}
	}*/
	return 0;
}

//Bマトリックスを求める関数
int Make_B_Matrix(int El_No, double B[D_MATRIX_SIZE][KIEL_SIZE], double Local_coord[DIMENSION],
	 double X[No_Control_point_ON_ELEMENT[Element_patch[El_No]]][DIMENSION], double *J, int Total_Control_Point){
	double a[DIMENSION][DIMENSION], b[DIMENSION][No_Control_point_ON_ELEMENT[Element_patch[El_No]]];


	int i,j,k;


		Jacobian( El_No, a, Local_coord, X, Total_Control_Point);


	*J = InverseMatrix_2D( a );
	//printf("B_Matri_J:%le\n",*J);
	if( *J <= 0 )return -999;

	for( i = 0; i < DIMENSION; i++ ){
		for( j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; j++ ){
			b[i][j] = 0.0;
			for( k = 0; k < DIMENSION; k++ ){
				b[i][j] += a[k][i] *dShape_func(j, k, Local_coord,El_No,Total_Control_Point);
				//b[i][j] += a[i][k] *dShape_func(j, k, Local_coord,El_No,Total_Control_Point);
			}
		}
	}
	for( i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++ ){
		B[0][2*i] = b[0][i];	B[0][2*i+1] = 0.0;
		B[1][2*i] = 0.0;		B[1][2*i+1] = b[1][i];
		B[2][2*i] = b[1][i];	B[2][2*i+1] = b[0][i];
	}
	/*for (i = 0; i < D_MATRIX_SIZE; i++)for (j = 0; j < KIEL_SIZE; j++) {
		printf("B[%d][%d]_B_mat:%le\n",i,j,B[i][j]);
	}*/
	return 0;
}
//応力歪マトリックス
int Make_D_Matrix_2D( double D[D_MATRIX_SIZE][D_MATRIX_SIZE], double E,double nu, int DM )
{
	int i,j;

	if(DM==0)//平面応力状態
	{
		//printf("E:%le nu:%le\n",E,nu);
		double Eone = E/(1.0-nu*nu);
		double D1[D_MATRIX_SIZE][D_MATRIX_SIZE] = {{ Eone, nu*Eone, 0},{nu*Eone, Eone, 0},{ 0, 0,(1-nu)/2*Eone} };

		for( i=0;i<D_MATRIX_SIZE;i++)
		for( j=0;j<D_MATRIX_SIZE;j++)
		D[i][j]=D1[i][j];
	}

	else if(DM==1)//平面ひずみ状態(2Dの場合はこっち)
	{
		//printf("E:%le nu:%le\n",E,nu);
		double Eone = E*(1.0-nu)/(1.0+nu)/(1.0-2*nu);
		double D1[D_MATRIX_SIZE][D_MATRIX_SIZE] = { {Eone, nu/(1.0-nu)*Eone,0},{ nu/(1.0-nu)*Eone, Eone, 0},{0, 0, (1-2*nu)/2/(1.0-nu)*Eone} };

		for( i=0;i<D_MATRIX_SIZE;i++)
		for( j=0;j<D_MATRIX_SIZE;j++)
		D[i][j]=D1[i][j];
	}

	else return ERROR;

	return 0;
}

//ガウスの数値積分法の中身
int BDBJ( double B[D_MATRIX_SIZE][KIEL_SIZE], double D[D_MATRIX_SIZE][D_MATRIX_SIZE], double J, double K_EL[KIEL_SIZE][KIEL_SIZE]){
	int i, j, k;
	double BD[KIEL_SIZE][D_MATRIX_SIZE];

	//[B]T[D][B]Jの計算
	for( i = 0; i < KIEL_SIZE; i++ ){
		for( j = 0; j < D_MATRIX_SIZE; j++ ){
			BD[i][j] = 0;
			for( k = 0; k < D_MATRIX_SIZE; k++ ){
				//printf("B[%d][%d]=%le D[%d][%d]=%le\n", k,i,B[k][i],k,j,D[k][j] );
				//printf("B[%d][%d]=%le\n", k,i,B[k][i]);
				BD[i][j] += B[k][i] * D[k][j];
				//printf("BD[%d][%d]=%e\n",i,j,BD[i][j] );
			}
		}
	}
	//for( j = 0; j < D_MATRIX_SIZE; j++ )for( i = 0; i < KIEL_SIZE; i++ )printf("B[%d][%d]=%le\n",j,i,B[j][i] );
	for( i = 0; i < KIEL_SIZE; i++ ){
		for( j = 0; j < KIEL_SIZE; j++ ){
			K_EL[i][j] = 0;
			for( k = 0; k < D_MATRIX_SIZE; k++ ){
				K_EL[i][j] +=  BD[i][k] *  B[k][j];
			}
			K_EL[i][j] *= J;
		}
	}
	return 0;
}

//要素合成マトリックス
int Make_K_EL(int El_No, double X[No_Control_point_ON_ELEMENT[Element_patch[El_No]]][DIMENSION], double K_EL[KIEL_SIZE][KIEL_SIZE], double E, double nu , int DM ,
	  int Total_Element, int Total_Control_Point){

	int i, j, k, l;

	double K1[KIEL_SIZE][KIEL_SIZE], B[D_MATRIX_SIZE][KIEL_SIZE];
	double D[D_MATRIX_SIZE][D_MATRIX_SIZE];
	double J = 0.0;
	double J_test = 0.0;
	double G = pow(0.6,0.5);
	double w[POW_Ng] = {(25.0/81.0),(40.0/81.0),(25.0/81.0),
						(40.0/81.0),(64.0/81.0),(40.0/81.0),
						(25.0/81.0),(40.0/81.0),(25.0/81.0)};
	double Gxi[POW_Ng][DIMENSION] = { {-G,-G},{0.0,-G},{G,-G},{-G,0.0},{0.0,0.0},{G,0.0},{-G,G},{0.0,G},{G,G} };
	//double Gxi[][POW_Ng][DIMENSION];


	for( i = 0; i < KIEL_SIZE; i++ ){
		for( j = 0; j < KIEL_SIZE; j++ ){
			K_EL[i][j] = 0.0;
		}
	}

	Make_D_Matrix_2D( D, E, nu ,DM);

	for( i = 0; i < POW_Ng; i++ ){

		//printf("i=%d\n",i );
		Make_B_Matrix( El_No, B,Gxi[i], X ,&J , Total_Control_Point);

		BDBJ( B, D, J, K1 );
		J_test += J;
		for( k = 0; k < KIEL_SIZE; k++ ){
			for( l = 0; l < KIEL_SIZE; l++ ){
				K_EL[k][l] += w[i] * K1[k][l];
			}
		}//printf("w[%d]=%f\n",i,w[i]);
	}printf("El_No:%d J_test=%e\n",El_No,J_test );
	//printf("G=%f\n",G );
	/*for ( k = 0; k < KIEL_SIZE; k++) {
		for ( l = 0; l < KIEL_SIZE; l++) {
			printf("K_EL[%d][%d]:%le\n",k,l,K_EL[k][l]);
		}
	}*/

	return 0;
}
///////////////////////////////////////////////////
//////////////歪と応力
//////////////////////////////////////////////////
void Make_Strain(double E, double nu, int Total_Element , int El_No, int Total_Control_Point){
	static double U[MAX_KIEL_SIZE];
	double B[D_MATRIX_SIZE][KIEL_SIZE],X[No_Control_point_ON_ELEMENT[Element_patch[El_No]]][DIMENSION],J;
	double G = pow(0.6,0.5);
	double Gxi[POW_Ng][DIMENSION] = { {-G,-G},{0.0,-G},{G,-G},{-G,0.0},{0.0,0.0},{G,0.0},{-G,G},{0.0,G},{G,G} };


	int N,e,i,j;
	//printf("Strain\n");
	for( e = 0; e < Total_Element; e++ ){
		//printf("\nElementNo.:%d\n",e);
		for( N = 0; N < POW_Ng; N++)
			for( i = 0; i < N_STRAIN; i ++ )
				Strain[e][N][i] = 0.0;
		//Bマトリックスと各要素の変位を取得
		//printf("El_No:%d,No_Control_point_ON_ELEMENT[%d]:%d\n", El_No,Element_patch[El_No],No_Control_point_ON_ELEMENT[Element_patch[El_No]]);
		for( i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[e]]; i++ ){
			for( j = 0; j < DIMENSION; j++ ){
				U[ i*DIMENSION +j ] = Displacement[ Controlpoint_of_Element[e][i]*DIMENSION + j ];
				X[i][j] = Node_Coordinate[ Controlpoint_of_Element[e][i] ][j];
			}
		}
		//歪
		for( N = 0; N < POW_Ng; N++ ){
			//printf("N:%d\n",N);
			Make_B_Matrix( e, B, Gxi[N], X ,&J , Total_Control_Point);
			for( i = 0; i < D_MATRIX_SIZE; i++ )
				for( j = 0; j < KIEL_SIZE; j++ ){
					 Strain[e][N][i] += B[i][j] * U[j] ;
					 //printf("B[%d][%d]_in_strain:%le * ",i,j,B[i][j]);
					 //if(e==1){
						 //printf("U[%d]=%le = %le\n",j,U[j],B[i][j]*U[j]);
					 //}
				 }
		}
	}
}
//応力
void Make_Stress_2D( double E, double nu ,int Total_Element, int DM){


	static double D[D_MATRIX_SIZE][D_MATRIX_SIZE];
	int e,i,j,k;
	Make_D_Matrix_2D( D, E, nu ,DM);


	for( e = 0; e < Total_Element; e++ ){
		for( k = 0; k < POW_Ng; k++)
			for( i = 0; i < N_STRESS; i ++ )
				Stress[e][k][i] = 0.0;
	for(k = 0; k < POW_Ng; k ++)
		for( i = 0; i < D_MATRIX_SIZE; i ++ )
			for( j = 0; j < D_MATRIX_SIZE; j++ )
				Stress[e][k][i] += D[i][j] * Strain[e][k][j];
	}
}

void Make_ReactionForce( int Total_Element,int Total_Control_Point, int El_No){
	int e, i, j, k,l,re;
	double B[D_MATRIX_SIZE][KIEL_SIZE],X[No_Control_point_ON_ELEMENT[Element_patch[El_No]]][DIMENSION],J;
	double w[POW_Ng] = {(25.0/81.0),(40.0/81.0),(25.0/81.0),
										  (40.0/81.0),(64.0/81.0),(40.0/81.0),
										  (25.0/81.0),(40.0/81.0),(25.0/81.0)};
	double G = pow(0.6,0.5);
	double Gxi[POW_Ng][DIMENSION] = { {-G,-G},{0.0,-G},{G,-G},{-G,0.0},{0.0,0.0},{G,0.0},{-G,G},{0.0,G},{G,G} };


	for( i = 0; i < Total_Control_Point* DIMENSION; i++ )
		ReactionForce[i] = 0.0;
		//printf("ReactionForce\n");
	for( re = 0; re < real_Total_Element; re++ ){
		e=real_element[re];
		//printf("%d,No_Control_point_ON_ELEMENT[%d]:%d\n",e,Element_patch[e], No_Control_point_ON_ELEMENT[Element_patch[e]]);
		//Bマトリックスを取得
		for( i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[e]]; i++ ){
			for( j = 0; j < DIMENSION; j++ ){
				X[i][j] = Node_Coordinate[ Controlpoint_of_Element[e][i] ][j];
				//printf("Node_Coordinate[%d][%d]:%le\n",Controlpoint_of_Element[e][i],j, Node_Coordinate[ Controlpoint_of_Element[e][i] ][j]);
				//printf("X[%d][%d]:%le\n",Controlpoint_of_Element[e][i],j,X[i][j] );
			}
		}
		for( k = 0; k < POW_Ng; k++ ){
			Make_B_Matrix( e, B, Gxi[k], X , &J, Total_Control_Point);
			for( j = 0; j < D_MATRIX_SIZE; j++ )
				for( l = 0; l < DIMENSION; l++)
					for( i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[e]]; i++ )
					ReactionForce[ Controlpoint_of_Element[e][i]*DIMENSION+l ]
							+= B[j][ i*DIMENSION+l ] * Stress[e][k][j] * w[k] * J;
							//printf("J:%le\n", J);
		}
	}
}
////////////////////////////////////////////////////////////////
//////////////////分布荷重//////////////////////////////////////
////////////////////////////////////////////////////////////////

//面番号ごとの節点番号の取り方の指定
/*void Force_Dis_NodeOfElement( int Number, int DistributeForce[MAX_N_DISTRIBUTE_FORCE][3], int ForceDis_NoE[No_Control_point_ON_ELEMENT[Element_patch[El_No]] ){
	int i;
	if( DistributeForce[Number][1] == 0 )
		for(i = 0; i < No_Control_point_ON_ELEMENT[iPatch]; i++ )
			ForceDis_NoE[i] = Controlpoint_of_Element[ DistributeForce[Number][0] ][i];
	else if( DistributeForce[Number][1] == 1){
		ForceDis_NoE[0] = Controlpoint_of_Element[ DistributeForce[Number][0] ][1];		ForceDis_NoE[1] = Controlpoint_of_Element[ DistributeForce[Number][0] ][2];
		ForceDis_NoE[2] = Controlpoint_of_Element[ DistributeForce[Number][0] ][3];		ForceDis_NoE[3] = Controlpoint_of_Element[ DistributeForce[Number][0] ][0];
	}
	else if( DistributeForce[Number][1] == 2){
		ForceDis_NoE[0] = Controlpoint_of_Element[ DistributeForce[Number][0] ][2];		ForceDis_NoE[1] = Controlpoint_of_Element[ DistributeForce[Number][0] ][3];
		ForceDis_NoE[2] = Controlpoint_of_Element[ DistributeForce[Number][0] ][0];		ForceDis_NoE[3] = Controlpoint_of_Element[ DistributeForce[Number][0] ][1];
	}
	else if( DistributeForce[Number][1] == 3){
		ForceDis_NoE[0] = Controlpoint_of_Element[ DistributeForce[Number][0] ][3];		ForceDis_NoE[1] = Controlpoint_of_Element[ DistributeForce[Number][0] ][0];
		ForceDis_NoE[2] = Controlpoint_of_Element[ DistributeForce[Number][0] ][1];		ForceDis_NoE[3] = Controlpoint_of_Element[ DistributeForce[Number][0] ][2];
	}

}
void Force_Dis( int Total_DistributeForce, int DistributeForce[MAX_N_DISTRIBUTE_FORCE][3], double Val_DistributeForce[MAX_N_DISTRIBUTE_FORCE],
					 int *Total_Load,int Load_Node_Dir[MAX_N_LOAD][2],double Value_of_Load[MAX_N_LOAD],int Total_Control_Point, int El_No, int *Total_Element ){


	int i,j, DF;
	static int ForceDistribute_Controlpoint_of_Element[MAX_NO_CCpoint_ON_ELEMENT];
	static double Out_Force_Distribute[MAX_NO_CCpoint_ON_ELEMENT];
	double a[DIMENSION][DIMENSION];
	static double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION];
	double G = pow( 0.6 , 0.5 );
	double J = 0.0,Sum=0.0;
	double w[DISTRIBUTE_FORCE_Ng] = {(5.0/9.0),(8.0/9.0),(5.0/9.0)};
	double Gxi[DISTRIBUTE_FORCE_Ng][DIMENSION] = { { 1.0, G },{1.0,0.0},{ 1.0, (-1.0)*G } };



	for( DF = 0; DF < Total_DistributeForce; DF++){
		//回転させた要素の節点座標の取得
		Force_Dis_NodeOfElement( DF, DistributeForce, ForceDistribute_Controlpoint_of_Element,El_No );

		for( i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++ ){
			for(j = 0; j < DIMENSION; j++ )
				X[i][j] = Node_Coordinate[ ForceDistribute_Controlpoint_of_Element[i] ][j];
			Out_Force_Distribute[i] = 0.0;
		}
		//各節点ごとの力の計算
		for( j = 0; j < DISTRIBUTE_FORCE_Ng; j++ ){
			Jacobian( a, Gxi[j], X , El_No,Total_Control_Point);
			J = pow(a[0][1],2)+pow(a[1][1],2);
			J = sqrt(J);
			for(i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++ )
				Out_Force_Distribute[i] +=Shape_func(i,Total_Control_Point,Gxi[j],El_No) * Val_DistributeForce[DF] * J * w[j];
		}
		//求めた値を荷重(Load)に追加
		for( i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++ ){
			if( Out_Force_Distribute[i] != 0 ){
				Load_Node_Dir[ *Total_Load ][0] = ForceDistribute_Controlpoint_of_Element[i];
				Load_Node_Dir[ *Total_Load ][1] = DistributeForce[DF][2];
				Value_of_Load[ *Total_Load ]  = Out_Force_Distribute[i];
				*Total_Load += 1;
				Sum += Out_Force_Distribute[i];
			}
		}
		printf("DistributeForce :Element= %d Factor= %d Direction= %d SumForce= %le\n",
			DistributeForce[DF][0],DistributeForce[DF][1],DistributeForce[DF][2],Sum);
		Sum=0.0;

	}

}*/

void Make_Parameter_z(int Total_Element, double E,double nu, int DM){
	int e,k;

	if(DM==0)
	{
	//Make_strain_z
	for(e=0;e<Total_Element;e++)
		for(k=0;k<POW_Ng;k++)
			Strain[e][k][3]=0.0;


	for(e=0;e<Total_Element;e++)
		for(k=0;k<POW_Ng;k++)
			Strain[e][k][3]=-1.0*nu/E*(Stress[e][k][0]+Stress[e][k][1]);
	}

	if(DM==1)
	{
	//Make_stree_z
	for(e=0;e<Total_Element;e++)
		for(k=0;k<POW_Ng;k++)
			Stress[e][k][3]=0.0;

	for(e=0;e<Total_Element;e++)
		for(k=0;k<POW_Ng;k++)
			Stress[e][k][3]=E*nu/(1.0+nu)/(1-2.0*nu)*(Strain[e][k][0]+Strain[e][k][1]);
	}

}

void element_coordinate(int Total_Element,int Total_Control_Point){
	int i,j,k,e,l;
	double element_edge[9][DIMENSION]={{-1.0,-1.0},{1.0,-1.0},{1.0,1.0},{-1.0,1.0},{0.0,-1.0},{1.0,0.0},{0.0,1.0},{-1.0,0.0},{0.0,0.0}};
	//double data_result_shape[2]={0.0};

	for ( e = 0; e < Total_Element; e++) {
		for ( k = 0; k < 9; k++) {
				double data_result_shape[2]={0.0};
			for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[e]]; i++) {
				for (j = 0; j < DIMENSION; j++) {
			data_result_shape[j] += Shape_func(i,Total_Control_Point,element_edge[k],e)*Node_Coordinate[Controlpoint_of_Element[e][i]][j];
				}
			}
			element_coordinate_Nopoint[l][0]=data_result_shape[0];
			element_coordinate_Nopoint[l][1]=data_result_shape[1];
			//printf("l;%d\n",l);
			/*printf("element_coordinate_Nopoint[%d][0]:%le element_coordinate_Nopoint[%d][1]:%le\n"
			,l,element_coordinate_Nopoint[l][0],l,element_coordinate_Nopoint[l][1]);*/
		l++;
			//printf("data_result_shape[0]:%le data_result_shape[1]:%le\n", data_result_shape[0],data_result_shape[1]);
		}
	}
	for (l = 0; l < 9*Total_Element; l++) {
		same_point_in_Element[l]=l;
	}
	for (l = 0; l < 9*Total_Element; l++) {
		for ( i = l-1; i >= 0 ; i--) {
		if (element_coordinate_Nopoint[l][0]==element_coordinate_Nopoint[i][0]&&element_coordinate_Nopoint[l][1]==element_coordinate_Nopoint[i][1]) {
			//printf("同じ座標の番号l:%d i:%d\n",l,i);
			same_point_in_Element[l]=i;
			}
		}
	}
}


void ourput_graph_2D(FILE *fp, int e, double element_gg, double element_ee,double data_result_shape_x, double data_result_shape_y, double data_result_disp_x, double data_result_disp_y){
	fp=fopen("NURBS/NURBS_disp.dat","a");
	fprintf(fp, "%d %20.13e  %20.13e %20.13e %20.13e %20.13e %20.13e\n",e,element_gg,element_ee,data_result_shape_x,data_result_shape_y,data_result_disp_x,data_result_disp_y);
	fclose(fp);
}

void Make_Strain_refine(double E, double nu, int Total_Element , int El_No, int Total_Control_Point,double Strain_ref[N_STRAIN],double element[DIMENSION]){
	static double U[MAX_KIEL_SIZE];
	double B[D_MATRIX_SIZE][KIEL_SIZE],X[No_Control_point_ON_ELEMENT[Element_patch[El_No]]][DIMENSION],J;



	int i,j;

	for( i = 0; i < N_STRAIN; i ++ )
	{
		Strain_ref[i] = 0.0;
	}
		//Bマトリックスと各要素の変位を取得
		for( i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++ ){
			for( j = 0; j < DIMENSION; j++ ){
				U[ i*DIMENSION +j ] = Displacement[ Controlpoint_of_Element[El_No][i]*DIMENSION + j ];
				X[i][j] = Node_Coordinate[ Controlpoint_of_Element[El_No][i] ][j];
			}
		}
		//歪

	//printf("N:%d\n",N);
	Make_B_Matrix( El_No, B, element, X ,&J , Total_Control_Point);
	for( i = 0; i < D_MATRIX_SIZE; i++ )
	{
		for( j = 0; j < KIEL_SIZE; j++ )
		{
			 Strain_ref[i] += B[i][j] * U[j];
		}
	}
}


void output_Stress_refine( int e, double element_gg, double element_ee,double data_result_shape[DIMENSION],double Stress_ref[N_STRAIN]){
	fp=fopen("Stress_ref/Stress_ref.dat","a");
	fprintf(fp, "%d %20.16lf  %20.16lf %20.16lf  %20.16lf %20.16lf %20.16lf %20.16lf %20.16lf\n",e,element_gg,element_ee,data_result_shape[0],data_result_shape[1],Stress_ref[0],Stress_ref[1],Stress_ref[2],Stress_ref[3]);
	fclose(fp);
}

void output_Stress_for_zarusoba(int p, double Stress_ref[N_STRESS]) {
	fp=fopen("colored_point/NURBS_stress_x.dat","a");
	fprintf(fp, "%d:%.16e\n",p,Stress_ref[0]);
	fclose(fp);

	fp=fopen("colored_point/NURBS_stress_y.dat","a");
	fprintf(fp, "%d:%.16e\n",p,Stress_ref[1]);
	fclose(fp);

	fp=fopen("colored_point/NURBS_stress_xy.dat","a");
	fprintf(fp, "%d:%.16e\n",p,Stress_ref[2]);
	fclose(fp);

	fp=fopen("colored_point/NURBS_stress_z.dat","a");
	fprintf(fp, "%d:%.16e\n",p,Stress_ref[3]);
	fclose(fp);

}
void output_disp_and_points_for_zarusoba(int p, double data_result_shape[DIMENSION],double data_result_disp[DIMENSION]) {
	//printf("stpr\n");
	//printf("p:%d %lf %lf %lf %lf\n", p,data_result_shape[0],data_result_shape[1],data_result_disp[0],data_result_disp[1]);
	if ((fp=fopen("colored_point/NURBS_points.txt","a"))!= NULL) {
		//printf("opencheck\n");
		//printf("%.16e %.16e\n",data_result_shape[0],data_result_shape[1]);
		fprintf(fp,"%.16e %.16e\n",data_result_shape[0],data_result_shape[1]);
		fclose(fp);
	}
	else{
		printf("cannot open NURBS_points.txt\n");
	}


	fp=fopen("colored_point/NURBS_disp_x.dat","a");
	fprintf(fp, "%d:%.16e\n",p,data_result_disp[0]);
	fclose(fp);

	fp=fopen("colored_point/NURBS_disp_y.dat","a");
	fprintf(fp, "%d:%.16e\n",p,data_result_disp[1]);
	fclose(fp);

		/*fp=fopen("colored_point/NURBS_disp_radius.dat","a");
		fprintf(fp, "%d:%.16e\n", p,pow((data_result_disp[0]*data_result_disp[0]+data_result_disp[1]*data_result_disp[1]),0.5));
		fclose(fp);*/
}

//積分点以上に細かく応力を計算するための関数

void Make_Stress_2D_refine( double E, double nu ,int Total_Element, int DM, int Total_Control_Point, double element[DIMENSION],double Stress_ref[N_STRESS],int El_No){


	static double D[D_MATRIX_SIZE][D_MATRIX_SIZE];
	static double Strain_ref[N_STRAIN];
	int i,j;

	Make_D_Matrix_2D( D, E, nu ,DM);

	Make_Strain_refine(E,nu,Total_Element,El_No,Total_Control_Point,Strain_ref,element);

	for( i = 0; i < N_STRESS; i ++ )
		Stress_ref[i] = 0.0;
	for( i = 0; i < D_MATRIX_SIZE; i ++ ){
		for( j = 0; j < D_MATRIX_SIZE; j++ ){
			Stress_ref[i] += D[i][j] * Strain_ref[j];
		}
	}
}


void calculate_Controlpoint_using_NURBS(double element[DIMENSION],int Total_Element,int Total_Control_Point,double E, double nu, int DM){
	int e,b,j,re,i;
	static double R;
    int p=0;
	//for(e=0; e < Total_Element; e++){



	for(re=0; re < real_Total_Element; re++){
		e=real_element[re];
		//printf("\n");
		//printf("Element_No:%d\n",e );
		double element_gg =0.0, element_ee=0.0, element_delta;

		int i_gg, i_ee/*, element_ndiv=10 */;

		No_points_for_colored_points=(element_ndiv+1)*(element_ndiv+1)*real_Total_Element;
		element_delta=2.0/element_ndiv;

		for(i_ee=0; i_ee < element_ndiv+1; i_ee++){
			for(i_gg=0; i_gg < element_ndiv+1; i_gg++){
				double data_result_shape[3]={0.0};
				double data_result_disp[3]={0.0};
				static double Stress_ref[N_STRESS];


				element_gg = -1.0 + element_delta*i_gg;
				element_ee = -1.0 + element_delta*i_ee;
				element[0]=element_gg;
				element[1]=element_ee;

		 //	printf("element_gg:%le\n",element_gg);
			 //printf("element_ee:%le\n",element_ee);


			 for (b = 0; b < No_Control_point_ON_ELEMENT[Element_patch[e]]; b++){
				 R=Shape_func(b,Total_Control_Point,element,e);
				 for (j = 0; j < DIMENSION; j++) {
					 //printf("%d %d %le\n",Controlpoint_of_Element[e][b],j,Displacement[Controlpoint_of_Element[e][b]*DIMENSION+j]);
					 //printf("%d %d %20.13le\n",Controlpoint_of_Element[e][b],j,Displacement[Controlpoint_of_Element[e][b]*DIMENSION+j]);
					 //printf("%d %d %le\n",Controlpoint_of_Element[e][b],j,Node_Coordinate[Controlpoint_of_Element[e][b]][j]);
					 data_result_disp[j] += R*Displacement[Controlpoint_of_Element[e][b]*DIMENSION+j];
					 data_result_shape[j] += R*Node_Coordinate[Controlpoint_of_Element[e][b]][j];
				 }
			 }

			 fp = fopen("shapefunc/shape_func_xi.dat", "a");
			 fprintf(fp,"%d %lf	%lf\t", e,element[0],Position_Data_param[0]);
			 for (i = 0; i < No_Control_point[Element_patch[e]][0]; i++) {
				 fprintf(fp, "%lf\t", Shape[0][i][Order[Element_patch[e]][0]]);
			 }
			 fprintf(fp,"\n");
			 fclose(fp);

			 fp = fopen("shapefunc/shape_func_eta.dat", "a");
			 fprintf(fp,"%d %lf	%lf\t", e,element[1],Position_Data_param[1]);
			 for (i = 0; i < No_Control_point[Element_patch[e]][1]; i++) {
				 fprintf(fp, "%lf\t", Shape[1][i][Order[Element_patch[e]][1]]);
			 }
			 fprintf(fp,"\n");
			 fclose(fp);
			 //data_result_shape_x[p]=data_result_shape[0];
			 //data_result_shape_y[p]=data_result_shape[1];
			 //p++;

			 //printf("\n");

			 ourput_graph_2D(fp, e,element_gg,element_ee,data_result_shape[0],data_result_shape[1],data_result_disp[0],data_result_disp[1]);
			 //data_result_shape_x[p]=data_result_shape[0];
			 //data_result_shape_y[p]=data_result_shape[1];
			 //data_result_disp_x[p]=data_result_disp[0];
			 //data_result_disp_y[p]=data_result_disp[1];

			 Make_Stress_2D_refine(E,nu,Total_Element,DM,Total_Control_Point,element,Stress_ref,e);
			 output_Stress_refine(e,element_gg,element_ee,data_result_shape,Stress_ref);

			 //output_Stress_for_zarusoba(p,Stress_ref);
			 //output_disp_and_points_for_zarusoba(p,data_result_shape,data_result_disp);
			 p++;

						//NURBS_points(fp,No_points_for_colored_points e,data_result_shape[0],data_result_shape[1])
					}
				}
			}
		}

void Gausspoint_coordinate(int Total_Element,int Total_Control_Point){
	int i,j,k,e;
	double G = pow(0.6,0.5);
	double Gxi[POW_Ng][DIMENSION] = { {-G,-G},{0.0,-G},{G,-G},{-G,0.0},{0.0,0.0},{G,0.0},{-G,G},{0.0,G},{G,G} };
	double R = 0.0;

	for ( e = 0; e < Total_Element; e++) {
		for ( k = 0; k < POW_Ng; k++) {
				double data_result_shape[2]={0.0};
			for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[e]]; i++) {
				R = Shape_func(i,Total_Control_Point,Gxi[k],e);
				for (j = 0; j < DIMENSION; j++) {
					//printf("i:%d Gxi[%d][%d]:%le\n", i,k,j,Gxi[k][i]);
			data_result_shape[j] += R *Node_Coordinate[Controlpoint_of_Element[e][i]][j];
				}
			}
			Gausspoint_coordinates[e][k][0]=data_result_shape[0];
			Gausspoint_coordinates[e][k][1]=data_result_shape[1];
		}
	}
}

/*void output_disp_and_points_for_extendmesh(int p, double data_result_shape[DIMENSION],double data_result_disp[DIMENSION]) {
	fp=fopen("new_zarusoba/extended_mesh.emsh","a");
	fprintf(fp, "%-.16lf %-.16lf",data_result_shape[0],data_result_shape[1]);
	fclose(fp);

	fp=fopen("new_zarusoba/displacement_x.dat","a");
	fprintf(fp, "%d:%-.16lf\n",p,data_result_disp[0]);
	fclose(fp);

	fp=fopen("new_zarusoba/displacement_y.dat","a");
	fprintf(fp, "%d:%-.16lf\n",p,data_result_disp[1]);
	fclose(fp);

		fp=fopen("colored_point/NURBS_disp_radius.dat","a");
		fprintf(fp, "%d:%-.16lf\n", p,pow((data_result_disp[0]*data_result_disp[0]+data_result_disp[1]*data_result_disp[1]),0.5));
		fclose(fp);
}*/


void calculate_extendmesh_using_NURBS(double element_emsh[DIMENSION],int Total_Element,int Total_Control_Point){
	int e,b,j,re;
	int p=Total_Control_Point;
	double R= 0.0;
	//for(e=0; e < Total_Element; e++){



	for(re=0; re < real_Total_Element; re++){
		e=real_element[re];
		//printf("\n");
		//printf("Element_No:%d\n",e );
		double element_gg =0.0, element_ee=0.0, element_delta_ee,element_delta_gg;

		int i_gg, i_ee,element_ndiv_ee = 5,element_ndiv_gg = 5;

		No_points_for_new_zarusoba=(element_ndiv_ee+1)*(element_ndiv_gg+1)*real_Total_Element;

		element_delta_ee=2.0/element_ndiv_ee;
		element_delta_gg=2.0/element_ndiv_gg;


		for(i_ee=0; i_ee < element_ndiv_ee+1; i_ee++){
			for(i_gg=0; i_gg < element_ndiv_gg+1; i_gg++){
				double data_result_shape[3]={0.0};
				double data_result_disp[3]={0.0};


				element_gg = -1.0 + element_delta_gg*i_gg;
				element_ee = -1.0 + element_delta_ee*i_ee;
				element_emsh[0]=element_gg;
				element_emsh[1]=element_ee;

		// printf("element_gg:%le\n",element_gg);
			//printf("element_ee:%le\n",element_ee);


						for (b = 0; b < No_Control_point_ON_ELEMENT[Element_patch[e]]; b++){
							R=Shape_func(b,Total_Control_Point,element_emsh,e);
							for (j = 0; j < DIMENSION; j++) {
							//printf("%d %d %le\n",Controlpoint_of_Element[e][b],j,Displacement[Controlpoint_of_Element[e][b]*DIMENSION+j]);
							//printf("%d %d %20.13le\n",Controlpoint_of_Element[e][b],j,Displacement[Controlpoint_of_Element[e][b]*DIMENSION+j]);
								data_result_disp[j] += R*Displacement[Controlpoint_of_Element[e][b]*DIMENSION+j];
								data_result_shape[j] += R*Node_Coordinate[Controlpoint_of_Element[e][b]][j];
					}
				}


								//data_result_shape_x[p]=data_result_shape[0];
								//data_result_shape_y[p]=data_result_shape[1];
								//p++;

				//printf("\n");

				//output_disp_and_points_for_extendmesh(p,data_result_shape,data_result_disp);

						data_result_shape_x_for_new_zarusoba[p]=data_result_shape[0];
						data_result_shape_y_for_new_zarusoba[p]=data_result_shape[1];
						data_result_disp_x_for_new_zarusoba[p]=data_result_disp[0];
						data_result_disp_y_for_new_zarusoba[p]=data_result_disp[1];
						p++;

						//NURBS_points(fp,No_points_for_colored_points e,data_result_shape[0],data_result_shape[1])

				}
			}
		}
	}

/*
////////////////////////////////////////////////////////////////
//////////////////AVS出力///////////////////////////////////////
////////////////////////////////////////////////////////////////

//節点座標と要素の節点番号の書き込み
void AVS_inputInp_Quad_4( int Total_Element, int Total_Control_Point ){
	int i;

	fprintf(fp,"%d	%d\n", Total_Control_Point, Total_Element);	//総節点数、総要素数
	for(i = 0; i < Total_Control_Point; i++ ){
		fprintf(fp,"%d	%e	%e	%e\n", i+1,Node_Coordinate[i][0],Node_Coordinate[i][1],0.0);
														//節点番号(1〜)、ＸＹＺ座標（２次元だとZ=0と記述
	}
	for( i = 0; i < Total_Element; i++ ){
		fprintf(fp,"%d	0	quad	", i+1);					//要素番号(1〜)、材料番号(0)、要素の形
		//要素タイプごとに順番を確認のとこ
		fprintf(fp,"%d	", Controlpoint_of_Element[i][2]+1);
		fprintf(fp,"%d	", Controlpoint_of_Element[i][3]+1);
		fprintf(fp,"%d	", Controlpoint_of_Element[i][0]+1);
		fprintf(fp,"%d	", Controlpoint_of_Element[i][1]+1);
		fprintf(fp,"\n");
	}

}
//計算結果のデータの書き込み
void AVS_inputAns_2D( int Total_Control_Point, int Total_Element ){

	int i,j,k, NodeDataCom = 6;
	int ElementDataCom = 12;	//N_STRAIN + N_STRESS + 4
	double Str;

	fprintf(fp,"%d	%d\n", NodeDataCom, ElementDataCom);	//各節点に存在するデータ数、各要素に存在するデータ数

	//節点のデータの書き込み（変位ＸＹＺ方向＋自由利用分XYZ）
	fprintf(fp,"%d	", NodeDataCom);						//節点データ成分数
	for( i = 0; i < NodeDataCom; i++ )
		fprintf(fp,"1	");								//各成分の構成数
	fprintf(fp,"\n");

	fprintf(fp,"DisX,\nDisY,\nDisZ,\nOpX,\nOpY,\nOpZ,\n");		//各節点データ成分のラベル
	for( i = 0; i < Total_Control_Point; i++ ){
		fprintf(fp,"%d	%e	%le	", i+1, Displacement[i*DIMENSION], Displacement[i*DIMENSION+1]);
		fprintf(fp,"%e	%le	%le	%le\n", 0.0, 0.0, 0.0, 0.0 );
	}

	//要素データの書き込み（応力＋自由利用２、歪＋自由利用２）
	fprintf(fp,"%d	", ElementDataCom);						//要素データ成分数
	for( i = 0; i < ElementDataCom; i++ )
		fprintf(fp,"1	");									//各成分の構成数
	fprintf(fp,"\n");

	fprintf(fp,"SigXX,\nSigYY,\nSigXY,\nSigZZ,\nSigOp1,\nSigOp2,\n");		//各節点データ成分のラベル
	fprintf(fp,"IpuXX,\nIpuYY,\nIpuXY,\nIpuZZ,\nIpuOp1,\nIpuOp2,\n");		//各節点データ成分のラベル

	for( i = 0; i < Total_Element; i++ ){
		fprintf(fp,"%d	", i+1);
		for( j = 0; j < N_STRESS; j++ ){
			Str = 0.0;
			for( k = 0; k < POW_Ng; k++ )		Str += Stress[i][k][j];
			fprintf(fp,"	%e",Str / (double)(POW_Ng) );
		}
		for( ; j < 6; j++ )				fprintf(fp,"	%e", 0.0 );

		for( j = 0; j < N_STRAIN; j++ ){
			Str = 0.0;
			for( k = 0; k < POW_Ng; k++ )
				Str += Strain[i][k][j];
			fprintf(fp,"	%e",Str / (double)(POW_Ng) );
		}
		for( ; j < 6; j++ )				fprintf(fp,"	%e", 0.0 );

		fprintf(fp,"\n");
	}
}

void Make_Output( int Total_Control_Point, int Total_Element ){
	int StepMax=1, StepNo=1;

	//AVS用のinpファイルの制作
	fp = fopen( "AVS/1_1_force_120_0122.inp", "w");
	fprintf(fp,"# AVS field file\n");					//注釈文（必ず先頭に「#」）
	fprintf(fp,"%d\n", StepMax);						//ステップ数の設定
	fprintf(fp,"data\n" );								//データの繰り返しタイプ	data,geom,data_geom

	fprintf(fp,"step%d\n", StepNo);					//ステップ番号
	AVS_inputInp_Quad_4( Total_Element, Total_Control_Point);	//節点座標と要素の取得
	AVS_inputAns_2D( Total_Control_Point, Total_Element );			//各種計算結果の取得
	printf("Finish Make_AVS_Step%d\n",StepNo);

	fclose(fp);
}
*/






/*
////////////////////////////////////////////////////////////////
//////////////////J積分/////////////////////////////////////////
////////////////////////////////////////////////////////////////



//B_xマトリックスを求める関数
int Make_B_x_Matrix_Quad_4(double B_x[DIMENSION][KIEL_SIZE], double Local_coord[DIMENSION], double X[No_Control_point_ON_ELEMENT][DIMENSION], double *J ){
	double a[DIMENSION][DIMENSION], b[DIMENSION][No_Control_point_ON_ELEMENT];
	int i,j,k;

	Jacobian_Quad_4( a, Local_coord, X );

	*J = InverseMatrix_2D( a );
	if( *J <= 0 )return -999;

	for( i = 0; i < DIMENSION; i++ ){
		for( j = 0; j < No_Control_point_ON_ELEMENT; j++ ){
			b[i][j] = 0.0;
			for( k = 0; k < DIMENSION; k++ ){
				b[i][j] += a[k][i] * dN_Quad_4( j, Local_coord, k);
			}
		}
	}
	for( i = 0; i < No_Control_point_ON_ELEMENT; i++ ){
		B_x[0][2*i] = b[0][i];	B_x[0][2*i+1] = 0.0;
		B_x[1][2*i] = 0.0;		B_x[1][2*i+1] = b[0][i];
	}
	return 0;
}







void Make_Strain_x_Quad_4(double E, double nu, int Total_Element){
	static double U[KIEL_SIZE];
	static double B_x[DIMENSION][KIEL_SIZE],X[No_Control_point_ON_ELEMENT][DIMENSION],J;
	double w[POW_Ng] = {1.0,1.0,1.0,1.0};
	double G = 1/pow(3,0.5);
	double Gxi[Total_Element][POW_Ng][DIMENSION] = { {{(1+G)/2,(1-G)/2},{(1+G)/2,(1+G)/2},{(1-G)/2,(1+G)/2},{(1-G)/2,(1-G)/2}},
	 {{(1+G)/2,(1-G)/2},{(1+G)/2,(1+G)/2},{(1-G)/2,(1+G)/2},{(1-G)/2,(1-G)/2}}};

	int N,e,i,j;

	for( e = 0; e < Total_Element; e++ ){
		for( N = 0; N < POW_Ng; N++)
			for( i = 0; i < DIMENSION; i ++ )
				Strain_x[e][N][i] = 0.0;
		//Bマトリックスと各要素の変位を取得
		for( i = 0; i < No_Control_point_ON_ELEMENT; i++ ){
			for( j = 0; j < DIMENSION; j++ ){
				U[ i*DIMENSION +j ] = Displacement[ Controlpoint_of_Element[e][i]*DIMENSION + j ];
				X[i][j] = Node_Coordinate[ Controlpoint_of_Element[e][i] ][j];
			}
		}
		//歪
		for( N = 0; N < POW_Ng; N++ ){
			Make_B_x_Matrix_Quad_4( B_x, Gxi[N], X ,&J );
			for( i = 0; i < DIMENSION; i++ )
				for( j = 0; j < KIEL_SIZE; j++ )
					 Strain_x[e][N][i] += B_x[i][j] * U[j] * w[N];
		}
	}
}












//エネルギーモーメンタムテンソルを求める関数
void Make_EMT(double E, double nu, int Total_Element){

	int i, j, k;
	double W;
	double W_x[DIMENSION];

	int K_D[DIMENSION] = {1, 0};
	double W_K_D[DIMENSION];
	double P_1j[DIMENSION];

	Make_Stress_2D(E, nu, Total_Element);
	Make_Strain_Quad_4(E, nu, Total_Element);
	Make_Strain_x_Quad_4(E, nu, Total_Element);

	W = 0.0;


	for( i = 0; i < DIMENSION; i++){
		W_x[i] = 0.0;
	}


	for( k = 0; k < POW_Ng; k++ ){
		for( i = 0; i < Total_Element; i++ ){
			for( j = 0; j < DIMENSION+1; j++ ){

				W += (1.0/2.0)*Stress[i][k][j]*Strain[i][k][j];
			}
		}
	}
	printf("\nW = %lf\n", W);


	for( i = 0; i < DIMENSION; i++ ){
		W_K_D[i] = W*K_D[i];

		printf("\nW_K_D[%d] = %lf\n", i, W_K_D[i]);
	}


	for( k = 0; k < POW_Ng; k++ ){
		for( i = 0; i < Total_Element; i++ ){
			for( j = 0; j < DIMENSION; j++ ){

				if(j == 0){
					W_x[j] += Stress[i][k][j]*Strain_x[i][k][j]+Stress[i][k][j+2]*Strain_x[i][k][j+1];
				}

				else if(j == 1){
					W_x[j] += Stress[i][k][j+1]*Strain_x[i][k][j-1]+Stress[i][k][j]*Strain_x[i][k][j];
				}
			}
		}
	}


	for( j = 0; j < DIMENSION; j++ ){
		printf("\nW_x[%d] = %le\n", j, W_x[j]);
	}


	for( i = 0; i < DIMENSION; i++ ){
		P_1j[i] = W_K_D[i] - W_x[i];

		printf("\nP_1j[%d] = %lf\n", i, P_1j[i]);
	}
}
*/

int  SerchForElement(int iPatch, int Total_Element, int iX, int iY){
	int iii;

	for(iii=0; iii<Total_Element; iii++){
		if (Element_patch[iii]==iPatch) {
			//(2019_10_10)printf("Check SerchForElement 1 iii = %d\n",iii);
			//(2019_10_10)printf("ENC[iii][0] = %d ENC[iii][1] = %d  iX = %d  iY = %d\n", ENC[iii][0], ENC[iii][1], iX, iY);
				if(iX == ENC[iii][0] && iY == ENC[iii][1]) goto loopend;
			/* iii --; */

			//(2019_10_10)printf("Check SerchForElement 2 iii = %d\n",iii);
			}
		}loopend:

	return (iii);
}

void Setting_Dist_Load_2D(int Total_Control_Point, int iPatch, int Total_Element, int iCoord, double val_Coord,
	 double Range_Coord[2], int type_load, double Coeff_Dist_Load[3]){
	int iii,jjj;
	int iDir_Element[MAX_N_KNOT],jDir_Element;
	int N_Seg_Load_Element_iDir=0, jCoord;
	int iRange_ele[2],jRange_ele;
	int iPos[2]={-10000, -10000}, jPos[2] = {-10000, -10000};
	double Coord_Seg_Load_Element_iDir[MAX_N_KNOT][2], Coord_Seg_Load_Element_jDir[2];
	int No_Element_for_Integration[MAX_N_KNOT], No_Element_For_Dist_Load;
	int iX, iY;
	//int Element_Integration;
	int iControlpoint[MAX_NO_CCpoint_ON_ELEMENT], ic, ig, NNG=3;
	double val_jCoord_Local;
	double GaussPt[3], Weight[3];
	double Gg = pow(3.0/5.0, 0.5);

	/* type_load: 0: Dist load in x direction
 * 	              1:              y direction
 * 	              2:              normal direciton */

	GaussPt[0] = - Gg; GaussPt[1] = 0.0; GaussPt[2] = Gg;
	Weight[0] = 5.0/9.0; Weight[1] = 8.0/9.0; Weight[2] = 5.0/9.0;

	/* iCoord=0: Load on Eta=Constant
	   iCoord=1: Load on Xi=Constant */
	if(iCoord == 0) jCoord = 1; if(iCoord == 1) jCoord = 0;

	/* val_Coord: Value of Eta or Xi of the line or surface to give the distributed load */


	/* Setting elements needed to computed the distributed load */

	for(iii=Order[iPatch][iCoord]; iii < No_knot[iPatch][iCoord]-Order[iPatch][iCoord]-1;iii++){
	double epsi=0.00000000001;
	/* iPos[0] = -10000; iPos[1] = -10000; jPos[0] = -10000; jPos[1] = -10000;*/
	//(2019_10_10)printf("Check1 iii = %d\n",iii);
	//(2019_10_10)printf("Check2 Position_Knots[iCoord][iii]= %f  Range_Coord[0] =%f Position_Knots[iCoord][iii+1] = %f\n",Position_Knots[iPatch][iCoord][iii],Range_Coord[0],Position_Knots[iPatch][iCoord][iii+1]);
	/*

		if(Position_Knots[iCoord][iii]-epsi <= Range_Coord[0] &&
			Position_Knots[iCoord][iii+1]+epsi > Range_Coord[0]) iPos[0] = iii;

		if(Position_Knots[iCoord][iii]-epsi <= Range_Coord[1] &&
                        Position_Knots[iCoord][iii+1]+epsi > Range_Coord[1]) iPos[1] = iii+1;
	*/
                if(Position_Knots[iPatch][iCoord][iii]-epsi <= Range_Coord[0]) iPos[0] = iii;
		if(Position_Knots[iPatch][iCoord][iii+1]-epsi <= Range_Coord[1]) iPos[1] = iii+1;

	}
	iRange_ele[0] = iPos[0] - Order[iPatch][iCoord]; iRange_ele[1] = iPos[1] - Order[iPatch][iCoord]-1;
	//(2019_10_10)printf("iPos[0] = %d iPos[1] = %d\n", iPos[0], iPos[1]);
	//(2019_10_10)printf("iRange_ele[0] = %d  iRange_ele[1] = %d\n",iRange_ele[0],iRange_ele[1]);

	if(iPos[0] < 0 || iPos[1] < 0) {
		printf("Error (Stop) iPos[0] = %d iPos[1] = %d\n", iPos[0], iPos[1]);
		exit(0);}

	for(jjj=Order[iPatch][jCoord]; jjj < No_knot[iPatch][jCoord]-Order[iPatch][jCoord]-1;jjj++){
	double epsi=0.00000000001;
	/* jjj=Order[jCoord]; */
		if(Position_Knots[iPatch][jCoord][jjj] - epsi <= val_Coord &&
			Position_Knots[iPatch][jCoord][jjj+1] + epsi> val_Coord) {jPos[0] = jjj; jPos[1] = jjj+1;
		 val_jCoord_Local = -1.0 + 2.0 * (val_Coord - Position_Knots[iPatch][jCoord][jjj]) /
			(Position_Knots[iPatch][jCoord][jjj+1] - Position_Knots[iPatch][jCoord][jjj]);
		}
		//(2019_06_13)printf("Check jjj count: jjj =  %d\n",jjj);
	}
	jRange_ele = jPos[0] - Order[iPatch][jCoord];
	//(2019_06_13)printf("jPos[0] = %d jPos[1] = %d  jRange_ele = %d val_jCoord_Local = %f\n",jPos[0],jPos[1],jRange_ele, val_jCoord_Local);

	 if(jPos[0] < 0 || jPos[1] < 0) {
                printf("Error (Stop) jPos[0] = %d jPos[1] = %d\n", jPos[0], jPos[1]);
                exit(0);}

	for(iii=iPos[0]; iii < iPos[1]; iii++){
		Coord_Seg_Load_Element_iDir[iii][0] = Position_Knots[iPatch][iCoord][iii+iPos[0]];
		Coord_Seg_Load_Element_iDir[iii][1] = Position_Knots[iPatch][iCoord][iii+iPos[0]+1];
		iDir_Element[N_Seg_Load_Element_iDir] = iii - Order[iPatch][iCoord];
		N_Seg_Load_Element_iDir++; }

		Coord_Seg_Load_Element_jDir[0] = Position_Knots[iPatch][iCoord][jPos[0]];
                Coord_Seg_Load_Element_jDir[1] = Position_Knots[iPatch][iCoord][jPos[1]];
                jDir_Element = jPos[0] - Order[iPatch][iCoord];
	iii=0;
	if(iCoord==1){
		int iX, iY;
		iX = jPos[0]-Order[iPatch][0];
		 //(2019_06_13)printf("Check iPos[0] = %d  iPos[1] = %d\n",iPos[0],iPos[1]);
		for(iY = iPos[0]-Order[iPatch][1]; iY < iPos[1]-Order[iPatch][1]; iY++){
			//(2019_06_13)printf("Check iY = %d\n",iY);
			No_Element_for_Integration[iii] = SerchForElement(iPatch, Total_Element, iX, iY);
			//(2019_06_13)printf("Check No_Element_for_Integration[%d] = %d\n",iii,No_Element_for_Integration[iii]);
			iii++;}
		}

	if(iCoord==0){
		int iX, iY;
                iY = jPos[0]-Order[iPatch][1];
		//(2019_06_13)printf("Check iPos[0] = %d  iPos[1] = %d\n",iPos[0],iPos[1]);
                for(iX = iPos[0]-Order[iPatch][0]; iX < iPos[1]-Order[iPatch][0]; iX++){
			//(2019_06_13)printf("Check iX = %d\n",iX);
                        No_Element_for_Integration[iii] = SerchForElement(iPatch, Total_Element, iX, iY);
			//(2019_06_13)printf("Check No_Element_for_Integration[%d] = %d\n",iii,No_Element_for_Integration[iii]);
                        iii++;}
		}
		No_Element_For_Dist_Load = iii;
	//(2019_06_13)printf("No_Element_For_Dist_Load = %d\n",No_Element_For_Dist_Load);


	/* Book keeping finished */

	for(iii=0; iii< No_Element_For_Dist_Load; iii++){  //B
		//(2019_06_13)printf("Check3 iii = %d\n",iii);
		//(2019_06_13)printf("Total_element_all_ID[No_Element_for_Integration[iii]] = %d\n No_Element_for_Integration[iii] = %d  iii = %d\n",
		//Total_element_all_ID[No_Element_for_Integration[iii]],No_Element_for_Integration[iii],iii);
		if(Total_element_all_ID[No_Element_for_Integration[iii]] == 1){  //A
			iX = ENC[No_Element_for_Integration[iii]][0];
			iY = ENC[No_Element_for_Integration[iii]][1];
			//(2019_06_13)printf("iX = %d  iY = %d\n",iX, iY);

			for(ic=0; ic < (Order[iPatch][0] + 1)*(Order[iPatch][1] + 1); ic++)
				iControlpoint[ic] =  Controlpoint_of_Element[No_Element_for_Integration[iii]][ic];

				for(ig = 0; ig < NNG; ig++){
					double Local_Coord[2],sfc,dxyzdge[3],detJ, XiEtaCoordParen, valDistLoad;
					int icc;
					Local_Coord[jCoord] = val_jCoord_Local;
					Local_Coord[iCoord] = GaussPt[ig];
					/*(2019_10_10)printf("ig = %d   Local_Coord[jCoord] = %f Local_Coord[iCoord] = %f\n"
						,ig,Local_Coord[jCoord],Local_Coord[iCoord]);*/

					ShapeFunc_from_paren(Local_Coord,iCoord,No_Element_for_Integration[iii]);
					XiEtaCoordParen = Position_Data_param[iCoord];
					/*(2019_10_10)printf("Check  Coeff_Dist_Load[0] = %f Coeff_Dist_Load[1] = %f  Coeff_Dist_Load[2] = %f  Position_Data_param[iCoord] = %f\n"
					, Coeff_Dist_Load[0], Coeff_Dist_Load[1], Coeff_Dist_Load[2],Position_Data_param[iCoord]);*/
					valDistLoad = Coeff_Dist_Load[0] + Coeff_Dist_Load[1] * XiEtaCoordParen
						+ Coeff_Dist_Load[2] * XiEtaCoordParen* XiEtaCoordParen;

					dxyzdge[0]=0.0; dxyzdge[1]=0.0;dxyzdge[2]=0.0;
					for(icc=0; icc <  (Order[iPatch][0] + 1)*(Order[iPatch][1] + 1); icc++){
					dxyzdge[0] += dShape_func(icc,iCoord,Local_Coord,No_Element_for_Integration[iii],
						Total_Control_Point) * Node_Coordinate[iControlpoint[icc]][0];
                                        dxyzdge[1] += dShape_func(icc,iCoord,Local_Coord,No_Element_for_Integration[iii],
                                                Total_Control_Point) * Node_Coordinate[iControlpoint[icc]][1];
					}

					detJ = sqrt(dxyzdge[0]*dxyzdge[0] + dxyzdge[1] * dxyzdge[1]);
					//(2019_06_13)printf("Check the value of detJ etc: detJ = %f dxyzdge[0] = %f dxyzdge[1] = %f\n",detJ,dxyzdge[0],dxyzdge[1]);
					if (type_load < 2) {
					for(ic = 0;  ic < (Order[iPatch][0] + 1)*(Order[iPatch][1] + 1); ic++){
						//printf("Order[%d][0];%d,Order[%d][1]:%d\n",iPatch,Order[iPatch][0],iPatch,Order[iPatch][1]);
						sfc=Shape_func(ic,Total_Control_Point, Local_Coord,
	                                                No_Element_for_Integration[iii]);
						Equivalent_Nodal_Force[iControlpoint[ic]][type_load] +=
							valDistLoad * sfc * detJ * Weight[ig];
					//(2019_06_13)printf("Check ic = %d sfc = %f   Weight[ig] = %f  valDistLoad = %f\n",ic,sfc,Weight[ig],valDistLoad);
					//(2019_06_13)printf("Equivalent_Nodal_Force[%d][%d]:%le\n",iControlpoint[ic],type_load,Equivalent_Nodal_Force[iControlpoint[ic]][type_load] );
						}
					}

					if (type_load == 2){
					double LoadDir[2];
					LoadDir[0] = dxyzdge[1]/detJ;
					LoadDir[1] = - dxyzdge[0]/detJ;
					for(ic = 0;  ic < (Order[iPatch][0] + 1)*(Order[iPatch][1] + 1); ic++){
                                                sfc=Shape_func(ic,Total_Control_Point, Local_Coord,
                                                        No_Element_for_Integration[iii]);
						Equivalent_Nodal_Force[iControlpoint[ic]][0] +=
                                                	LoadDir[0] * valDistLoad * sfc * detJ * Weight[ig];
						Equivalent_Nodal_Force[iControlpoint[ic]][1] +=
                                                        LoadDir[1] * valDistLoad * sfc * detJ * Weight[ig];
						}
					}
				}
			} //A
		} //B
}
