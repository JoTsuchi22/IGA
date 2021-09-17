/***************************
 * s-IGA用のNURBS_input
 * ローカルパッチについて、グローバルと重ね合わせた結果を出力
 * 問題点（グローバルパッチのパッチ表面において正確なひずみ算出が困難）の解決を目指す:済
 * グローバルの要素境界における結合結果算出は可能か要確認
 * 
 *************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include "constant.h"
#include "NURBS_calc.h"

#define DBL_MAX          1.7976931348623158e+308 // max value

int stress_type_flag = 1;			//0:平面応力 1:平面ひずみ
double E;							//ヤング率(GPa)
double nu;							//ポアソン比(-)
int patch_n;						//パッチ数
int cntl_p_n;						//コントロールポイント数
int order_xi[MAX_PATCHES];		//ξ基底関数の次数(p)
int order_eta[MAX_PATCHES];		//η基底関数の次数(p)
int knot_n_xi[MAX_PATCHES];		//ξノットベクトルの数(n+p+1)
int knot_n_eta[MAX_PATCHES];		//ηノットベクトルの数(n+p+1)
int cntl_p_n_xi[MAX_PATCHES];	//ξ方向コントロールポイント数(n)
int cntl_p_n_eta[MAX_PATCHES];	//η方向コントロールポイント数(n)

double knot_vec_xi[MAX_PATCHES][MAX_KNOTS];		//ξノットベクトル
double knot_vec_eta[MAX_PATCHES][MAX_KNOTS];	//ηノットベクトル
double cntl_px[MAX_PATCHES][MAX_CNRL_P];		//コントロールポイントx座標
double cntl_py[MAX_PATCHES][MAX_CNRL_P];		//コントロールポイントy座標
double disp_cntl_px[MAX_PATCHES][MAX_CNRL_P];	//コントロールポイント上のx方向変位
double disp_cntl_py[MAX_PATCHES][MAX_CNRL_P];	//コントロールポイント上のy方向変位
double weight[MAX_PATCHES][MAX_CNRL_P];			//重み

double coord_x[MAX_POINTS][MAX_POINTS];		//メッシュx座標
double coord_y[MAX_POINTS][MAX_POINTS];		//メッシュy座標
double dxi_x[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];		// ∂x/∂ξ
double dxi_y[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];		// ∂y/∂ξ
double deta_x[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];		// ∂x/∂η
double deta_y[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];		// ∂y/∂η

double disp_x[MAX_POINTS][MAX_POINTS];			//x方向変位
double disp_y[MAX_POINTS][MAX_POINTS];			//y方向変位
double dxi_disp_x[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	// ∂u/∂ξ
double dxi_disp_y[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	// ∂v/∂ξ
double deta_disp_x[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	// ∂u/∂η
double deta_disp_y[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	// ∂v/∂η

double strain_xx[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//x方向ひずみ
double strain_yy[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//y方向ひずみ
double strain_xy[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//剪断ひずみ

double stress_xx[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//x方向垂直応力
double stress_yy[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//y方向垂直応力
double stress_xy[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//剪断応力

int fields_flag = 1;		//s-IGAのためのNURBS_inputでは変位データは必ず読み込ませる
int division_ele_xi;		//ξ方向の一要素あたりの分割数
int division_ele_eta;		//η方向の一要素あたりの分割数
int division_n_xi;		//ξ方向の表示する点の数
int division_n_eta;		//η方向の表示する点の数
int element_n_xi;			//ξ方向要素数
int element_n_eta;		//η方向要素数

int temp_index[MAX_PATCHES][MAX_CNRL_P];
double temp_cntl_px[MAX_CNRL_P];
double temp_cntl_py[MAX_CNRL_P];
double temp_weight[MAX_CNRL_P];
double temp_disp_x[MAX_CNRL_P];
double temp_disp_y[MAX_CNRL_P];

//for s-IGA
int n_patch_glo;	//グローバルメッシュ上のパッチ数
int n_patch_loc;	//ローカルメッシュ上のパッチ数
int glo_cntl_p_n;	//グローバルメッシュ上のコントロールポイント数
int loc_cntl_p_n;	//ローカルメッシュ上のコントロールポイント数

//for graph
int graph_patch_n;	//グラフ作成用出力ファイル内のパッチ番号

FILE *fp;	// FILE型構造体

void GetLocData (char *filename_input_loc)
{
	double temp;

	//必要なのはローカルのパッチ数とコントロールポイント数
	printf("Start Get Local Data\n\n");
	fp = fopen(filename_input_loc, "r");

	fscanf(fp, "%lf%lf", &temp, &temp);
	fscanf(fp, "\n");
	fscanf(fp, "%d", &n_patch_loc);
	fscanf(fp, "\n");
	fscanf(fp, "%d", &loc_cntl_p_n);
	fscanf(fp, "\n");
	fclose(fp);
	printf("patches(in local):% d\n",n_patch_loc);
	printf("control points(in local):% d\n",loc_cntl_p_n);
	printf("\nFinish Get Local Data\n\n");

}

void ReadFile (char *filename_input, char *filename_disp) {
	int i, j;
	double temp1, temp2, temp3;
	int temp_int;

	printf("Start Reading input\n\n");
	fp = fopen(filename_input, "r");

	fscanf(fp, "%lf%lf", &E, &nu);
	fscanf(fp, "\n");
	printf("E nu: % 1.4e % 1.4e\n", E, nu);

	fscanf(fp, "%d", &patch_n);
	fscanf(fp, "\n");
	printf("patches: %d \n", patch_n);
	if (patch_n > MAX_PATCHES) {
		printf("Error!!\n");
		printf("Too many patches!\n"
			   "Maximum of patches is %d (Now %d)\n"
			   "\n", MAX_PATCHES, patch_n);
		exit(1);
	}

	fscanf(fp, "%d", &cntl_p_n);
	fscanf(fp, "\n");
	printf("total control points:%d \n", cntl_p_n);
	if (cntl_p_n > MAX_CNRL_P) {
		printf("Error!!\n");
		printf("Too many control points!\n"
			   "Maximum of control points is %d (Now %d)\n"
			   "\n", MAX_CNRL_P, cntl_p_n);
		exit(1);
	}

	for (i = 0; i < patch_n; i++) {
		fscanf(fp, "%d%d", &order_xi[i], &order_eta[i]);
		fscanf(fp, "\n");
		printf("order %d: %d %d\n", i, order_xi[i], order_eta[i]);
		if (order_xi[i] > MAX_ORDER) {
			printf("Error!!\n");
			printf("Order too big at xi!\n"
				   "Maximum of order is %d (Now %d at patch %d)\n"
				   "\n", MAX_ORDER, order_xi[i], i);
			exit(1);
		}
		if (order_eta[i] > MAX_ORDER) {
			printf("Error!!\n");
			printf("Order too big at eta!\n"
				   "Maximum of order is %d (Now %d at patch %d)\n"
				   "\n", MAX_ORDER, order_eta[i], i);
			exit(1);
		}
	}

	for (i = 0; i < patch_n; i++) {
		fscanf(fp, "%d%d", &knot_n_xi[i], &knot_n_eta[i]);
		fscanf(fp, "\n");
		printf("knots %d: %d %d\n", i, knot_n_xi[i], knot_n_eta[i]);
		if (knot_n_xi[i] > MAX_KNOTS) {
			printf("Error!!\n");
			printf("Knot vector too long at xi!\n"
				   "Maximum of knot vector is %d (Now %d at patch %d)\n"
				   "\n", MAX_KNOTS, knot_n_xi[i], i);
			exit(1);
		}
		if (knot_n_eta[i] > MAX_KNOTS) {
			printf("Error!!\n");
			printf("Knot vector too long at eta!\n"
				   "Maximum of knot vector is %d (Now %d at patch %d)\n"
				   "\n", MAX_KNOTS, knot_n_eta[i], i);
			exit(1);
		}
	}

	for (i = 0; i < patch_n; i++) {
		fscanf(fp, "%d%d", &cntl_p_n_xi[i], &cntl_p_n_eta[i]);
		printf("control points %d: %d %d\n",
			   i, cntl_p_n_xi[i], cntl_p_n_eta[i]);
		fscanf(fp, "\n");
	}
	printf("\n");


	for (i = 0; i < patch_n; i++) {
		for (j = 0; j < cntl_p_n_xi[i] * cntl_p_n_eta[i]; j++) {
			fscanf(fp, "%d", &temp_index[i][j]);
			printf("%d ", temp_index[i][j]);
		}
		fscanf(fp, "\n");
		printf("\n");
	}
	printf("\n");

	fscanf(fp, "%lf%lf%lf", &temp1, &temp2, &temp3);
	fscanf(fp, "\n");

	for (i = 0; i < patch_n; i++) {
		for (j = 0; j < knot_n_xi[i]; j++) {
			fscanf(fp, "%le", &knot_vec_xi[i][j]);
			printf("%f\t", knot_vec_xi[i][j]);
		}
		printf("\n");
		for (j = 0; j < knot_n_eta[i]; j++) {
			fscanf(fp, "%le", &knot_vec_eta[i][j]);
			printf("%f\t", knot_vec_eta[i][j]);
		}
		printf("\n");
	}
	printf("\n");

	for (i = 0; i < cntl_p_n; i++) {
		fscanf(fp, "%d%le%le%le",
			   &temp_int,
			   &temp_cntl_px[i], &temp_cntl_py[i], &temp_weight[i]);
		printf("%d\t%f\t%f\t%f\n",
			   temp_int,
			   temp_cntl_px[i], temp_cntl_py[i], temp_weight[i]);
	}
	printf("\n");

	for (i = 0; i < patch_n; i++) {
		for (j = 0; j < cntl_p_n_xi[i] * cntl_p_n_eta[i]; j++) {
			cntl_px[i][j] = temp_cntl_px[temp_index[i][j]];
			cntl_py[i][j] = temp_cntl_py[temp_index[i][j]];
			weight[i][j] = temp_weight[temp_index[i][j]];
			printf("%d\t%f\t%f\t%f\n",
				   temp_index[i][j],
				   cntl_px[i][j], cntl_py[i][j], weight[i][j]);
		}
		printf("\n");
	}
	fclose(fp);
	printf("End Reading input\n\n");

	if (fields_flag) {
		printf("Start Reading displacement\n\n");
		fp = fopen(filename_disp, "r");
		char buff[256];

		fscanf(fp, "%s", buff);
		fscanf(fp, "%s", buff);

		for (i = 0; i < cntl_p_n; i++) {
			fscanf(fp, "%d:%le%le",
				   &temp_int, &temp_disp_x[i], &temp_disp_y[i]);
			printf("%d\t%1.6e\t%1.6e\n",
				   temp_int, temp_disp_x[i], temp_disp_y[i]);
		}
		printf("\n");

		for (i = 0; i < patch_n; i++) {
			for (j = 0; j < cntl_p_n_xi[i] * cntl_p_n_eta[i]; j++) {
				disp_cntl_px[i][j] = temp_disp_x[temp_index[i][j]];
				disp_cntl_py[i][j] = temp_disp_y[temp_index[i][j]];
				printf("%d\t%f\t%f\t%f\n",
					   temp_index[i][j], cntl_px[i][j], cntl_py[i][j], weight[i][j]);
			}
			printf("\n");
		}
		fclose(fp); // ファイルを閉じる
		printf("End Reading displpacement\n\n");
	}

	fp = fopen("Displacement_loc.dat", "w");
	glo_cntl_p_n = cntl_p_n - loc_cntl_p_n;
	fprintf(fp, "label=Displacement\n"
				"num_items=%d\n\n", loc_cntl_p_n);
	for (i = 0; i < loc_cntl_p_n; i++)
	{
		fprintf(fp, "%d:	%le %le \n", 
				i, temp_disp_x[i + glo_cntl_p_n], temp_disp_y[i + glo_cntl_p_n]);
	}
	fclose(fp);
}

int CalcXiEtaByNR(double px, double py,
                  double *knot_vec_xi, double *knot_vec_eta,
                  double *cntl_px, double *cntl_py,
                  double *disp_cntl_px, double *disp_cntl_py,
                  int cntl_p_n_xi, int cntl_p_n_eta,
                  double *weight, int order_xi, int order_eta,
                  double *output_xi, double *output_eta,
				  double *disp_x_glo, double *disp_y_glo,
                  double *strain_xx_glo, double *strain_yy_glo, double *strain_xy_glo) {
	double temp_xi, temp_eta;
	double temp_x, temp_y;
	double temp_matrix[2][2];
	double temp_dxi, temp_deta;
	double temp_tol_x = DBL_MAX;
	double temp_tol_y = DBL_MAX;

	(*output_xi) = 0;
	(*output_eta) = 0;

	int i;
	double tol = 10e-8;

	temp_xi = knot_vec_xi[0] + knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = knot_vec_eta[0] + knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;
	//printf("% 1.8e % 1.8e\n", temp_xi, temp_eta);
	
	for (i = 0; i < 1000; i++) {
		rNURBS_surface(knot_vec_xi, knot_vec_eta,
		               cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
		               weight, order_xi, order_eta,
		               temp_xi, temp_eta,
		               &temp_x, &temp_y,
		               &temp_matrix[0][0], &temp_matrix[0][1],
		               &temp_matrix[1][0], &temp_matrix[1][1]);

		temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;
		
		//収束した場合////////////////////////////////////////////////////////////////
		//if (temp_tol_x < tol && temp_tol_y < tol) {
        if (temp_tol_x + temp_tol_y < tol) {
			printf("rNURBS\n");
			if (temp_xi == knot_vec_xi[0] || temp_eta == knot_vec_eta[0])
			{
				break;
			}
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;

			int knot_n_xi = cntl_p_n_xi + order_xi + 1;
			int knot_n_eta = cntl_p_n_eta + order_eta + 1;
			double dtilda_xi = 0.0;
			double dtilda_eta = 0.0;
			double disp_x, disp_y;
			double dxi_x, deta_x, dxi_y, deta_y;
			double dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y;
			double temp_matrix2[2][2];
			double temp;
			double strain_xx, strain_yy, strain_xy;
			double stress_xx, stress_yy, stress_xy;

			for (i = 0; i < knot_n_xi; i++) {
				if ( knot_vec_xi[i] < temp_xi && temp_xi <= knot_vec_xi[i + 1]) {
					dtilda_xi = ( knot_vec_xi[i + 1] - knot_vec_xi[i] ) / 2.0;
					printf("xi%f\n", dtilda_xi);
					break;
				}
			}
			for (i = 0; i < knot_n_eta; i++) {
				if ( knot_vec_eta[i] < temp_eta && temp_eta <= knot_vec_eta[i + 1]) {
					dtilda_eta = ( knot_vec_eta[i + 1] - knot_vec_eta[i] ) / 2.0;
					printf("eta%f\n", dtilda_eta);
					break;
				}
			}

			rNURBS_surface(knot_vec_xi, knot_vec_eta,
			               cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
			               weight, order_xi, order_eta,
			               temp_xi, temp_eta,
			               &temp, &temp,
			               &dxi_x, &deta_x,
			               &dxi_y, &deta_y);
			printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
				   dxi_x, deta_x, dxi_y, deta_y);

			rNURBS_surface(knot_vec_xi, knot_vec_eta,
			               disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
			               weight, order_xi, order_eta,
			               temp_xi, temp_eta,
			               &disp_x, &disp_y,
			               &dxi_disp_x, &deta_disp_x,
			               &dxi_disp_y, &deta_disp_y);
			printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
				   dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y);

			temp_matrix2[0][0] = dxi_x * dtilda_xi;
			temp_matrix2[0][1] = dxi_y * dtilda_xi;
			temp_matrix2[1][0] = deta_x * dtilda_eta;
			temp_matrix2[1][1] = deta_y * dtilda_eta;

			InverseMatrix_2D(temp_matrix2);

			strain_xx = temp_matrix2[0][0] * dtilda_xi
			            * dxi_disp_x
			            + temp_matrix2[0][1] * dtilda_eta
			            * deta_disp_x;
			strain_yy = temp_matrix2[1][0] * dtilda_xi
			            * dxi_disp_y
			            + temp_matrix2[1][1] * dtilda_eta
			            * deta_disp_y;
			strain_xy = temp_matrix2[1][0] * dtilda_xi
			            * dxi_disp_x
			            + temp_matrix2[1][1] * dtilda_eta
			            * deta_disp_x
			            + temp_matrix2[0][0] * dtilda_xi
			            * dxi_disp_y
			            + temp_matrix2[0][1] * dtilda_eta
			            * deta_disp_y;

			double D_matrix[3][3] = {{0.0}};
			int stress_type_flag = 1;
			if (stress_type_flag == 0) { //平面応力状態
				temp = E * (1.0 - nu * nu);
				D_matrix[0][0] = temp;
				D_matrix[0][1] = nu * temp;
				D_matrix[1][0] = nu * temp;
				D_matrix[1][1] = temp;
				D_matrix[2][2] = (1.0 - nu) / 2.0 * temp;
			} else if (stress_type_flag == 1) { //平面ひずみ状態(2Dの場合はこっち)
				temp = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
				D_matrix[0][0] = temp;
				D_matrix[0][1] = nu / (1.0 - nu) * temp;
				D_matrix[1][0] = nu / (1.0 - nu) * temp;
				D_matrix[1][1] = temp;
				D_matrix[2][2] = (1.0 - 2.0 * nu) / 2.0 / (1.0 - nu) * temp;
			}

			stress_xx = D_matrix[0][0] * strain_xx
			            + D_matrix[0][1] * strain_yy;
			stress_yy = D_matrix[1][0] * strain_xx
			            + D_matrix[1][1] * strain_yy;
			stress_xy = D_matrix[2][2] * strain_xy;

			temp = sqrt(disp_x * disp_x + disp_y * disp_y);
			printf("x:   % 1.8e\n", px);
			printf("y:   % 1.8e\n", py);
			printf("xi:  % 1.8e\n", temp_xi);
			printf("eta: % 1.8e\n", temp_eta);
			printf("Displacement x: % 1.8e\n", disp_x);
			printf("Displacement y: % 1.8e\n", disp_y);
			printf("Displacement  : % 1.8e\n", temp);
			printf("Strain xx: % 1.8e\n", strain_xx);
			printf("Strain yy: % 1.8e\n", strain_yy);
			printf("Strain xy: % 1.8e\n", strain_xy);
			printf("Stress xx: % 1.8e\n", stress_xx);
			printf("Stress yy: % 1.8e\n", stress_yy);
			printf("Stress xy: % 1.8e\n", stress_xy);

			(*disp_x_glo) = disp_x;
			(*disp_y_glo) = disp_y;

			(*strain_xx_glo) = strain_xx;
			(*strain_yy_glo) = strain_yy;
			(*strain_xy_glo) = strain_xy;

			return i;
		}

		InverseMatrix_2D(temp_matrix);

		temp_dxi = temp_matrix[0][0] * (px - temp_x)
		           + temp_matrix[0][1] * (py - temp_y);
		temp_deta = temp_matrix[1][0] * (px - temp_x)
		            + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_deta;
		if (temp_xi < knot_vec_xi[0])
			temp_xi = knot_vec_xi[0];
		if (temp_xi > knot_vec_xi[cntl_p_n_xi + order_xi])
			temp_xi = knot_vec_xi[cntl_p_n_xi + order_xi];
		if (temp_eta < knot_vec_eta[0])
			temp_eta = knot_vec_eta[0];
		if (temp_eta > knot_vec_eta[cntl_p_n_eta + order_eta])
			temp_eta = knot_vec_eta[cntl_p_n_eta + order_eta];

		//temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
		//printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
	}
	temp_xi = knot_vec_xi[0] + knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = knot_vec_eta[0] + knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;
	for (i = 0; i < 1000; i++) {
		lNURBS_surface(knot_vec_xi, knot_vec_eta,
		               cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
		               weight, order_xi, order_eta,
		               temp_xi, temp_eta,
		               &temp_x, &temp_y,
		               &temp_matrix[0][0], &temp_matrix[0][1],
		               &temp_matrix[1][0], &temp_matrix[1][1]);

		temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;
		
		//収束した場合////////////////////////////////////////////////////////////////
		//if (temp_tol_x < tol && temp_tol_y < tol) {
        if (temp_tol_x + temp_tol_y < tol) {
			printf("lNURBS\n");		
			if (temp_xi == knot_vec_xi[cntl_p_n_xi + order_xi] || temp_eta == knot_vec_eta[cntl_p_n_eta + order_eta])
			{
				break;
			}
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;

			int knot_n_xi = cntl_p_n_xi + order_xi + 1;
			int knot_n_eta = cntl_p_n_eta + order_eta + 1;
			double dtilda_xi = 0.0;
			double dtilda_eta = 0.0;
			double disp_x, disp_y;
			double dxi_x, deta_x, dxi_y, deta_y;
			double dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y;
			double temp_matrix2[2][2];
			double temp;
			double strain_xx, strain_yy, strain_xy;
			double stress_xx, stress_yy, stress_xy;

			for (i = 0; i < knot_n_xi; i++) {
				if ( knot_vec_xi[i] <= temp_xi && temp_xi < knot_vec_xi[i + 1]) {
					dtilda_xi = ( knot_vec_xi[i + 1] - knot_vec_xi[i] ) / 2.0;
					printf("%f\n", dtilda_xi);
					break;
				}
			}
			for (i = 0; i < knot_n_eta; i++) {
				if ( knot_vec_eta[i] <= temp_eta && temp_eta < knot_vec_eta[i + 1]) {
					dtilda_eta = ( knot_vec_eta[i + 1] - knot_vec_eta[i] ) / 2.0;
					printf("%f\n", dtilda_eta);
					break;
				}
			}

			lNURBS_surface(knot_vec_xi, knot_vec_eta,
			               cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
			               weight, order_xi, order_eta,
			               temp_xi, temp_eta,
			               &temp, &temp,
			               &dxi_x, &deta_x,
			               &dxi_y, &deta_y);
			printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
				   dxi_x, deta_x, dxi_y, deta_y);

			lNURBS_surface(knot_vec_xi, knot_vec_eta,
			               disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
			               weight, order_xi, order_eta,
			               temp_xi, temp_eta,
			               &disp_x, &disp_y,
			               &dxi_disp_x, &deta_disp_x,
			               &dxi_disp_y, &deta_disp_y);
			printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
				   dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y);

			temp_matrix2[0][0] = dxi_x * dtilda_xi;
			temp_matrix2[0][1] = dxi_y * dtilda_xi;
			temp_matrix2[1][0] = deta_x * dtilda_eta;
			temp_matrix2[1][1] = deta_y * dtilda_eta;

			InverseMatrix_2D(temp_matrix2);

			strain_xx = temp_matrix2[0][0] * dtilda_xi
			            * dxi_disp_x
			            + temp_matrix2[0][1] * dtilda_eta
			            * deta_disp_x;
			strain_yy = temp_matrix2[1][0] * dtilda_xi
			            * dxi_disp_y
			            + temp_matrix2[1][1] * dtilda_eta
			            * deta_disp_y;
			strain_xy = temp_matrix2[1][0] * dtilda_xi
			            * dxi_disp_x
			            + temp_matrix2[1][1] * dtilda_eta
			            * deta_disp_x
			            + temp_matrix2[0][0] * dtilda_xi
			            * dxi_disp_y
			            + temp_matrix2[0][1] * dtilda_eta
			            * deta_disp_y;

			double D_matrix[3][3] = {{0.0}};
			int stress_type_flag = 1;
			if (stress_type_flag == 0) { //平面応力状態
				temp = E * (1.0 - nu * nu);
				D_matrix[0][0] = temp;
				D_matrix[0][1] = nu * temp;
				D_matrix[1][0] = nu * temp;
				D_matrix[1][1] = temp;
				D_matrix[2][2] = (1.0 - nu) / 2.0 * temp;
			} else if (stress_type_flag == 1) { //平面ひずみ状態(2Dの場合はこっち)
				temp = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
				D_matrix[0][0] = temp;
				D_matrix[0][1] = nu / (1.0 - nu) * temp;
				D_matrix[1][0] = nu / (1.0 - nu) * temp;
				D_matrix[1][1] = temp;
				D_matrix[2][2] = (1.0 - 2.0 * nu) / 2.0 / (1.0 - nu) * temp;
			}

			stress_xx = D_matrix[0][0] * strain_xx
			            + D_matrix[0][1] * strain_yy;
			stress_yy = D_matrix[1][0] * strain_xx
			            + D_matrix[1][1] * strain_yy;
			stress_xy = D_matrix[2][2] * strain_xy;

			temp = sqrt(disp_x * disp_x + disp_y * disp_y);
			printf("x:   % 1.8e\n", px);
			printf("y:   % 1.8e\n", py);
			printf("xi:  % 1.8e\n", temp_xi);
			printf("eta: % 1.8e\n", temp_eta);
			printf("Displacement x: % 1.8e\n", disp_x);
			printf("Displacement y: % 1.8e\n", disp_y);
			printf("Displacement  : % 1.8e\n", temp);
			printf("Strain xx: % 1.8e\n", strain_xx);
			printf("Strain yy: % 1.8e\n", strain_yy);
			printf("Strain xy: % 1.8e\n", strain_xy);
			printf("Stress xx: % 1.8e\n", stress_xx);
			printf("Stress yy: % 1.8e\n", stress_yy);
			printf("Stress xy: % 1.8e\n", stress_xy);

			(*disp_x_glo) = disp_x;
			(*disp_y_glo) = disp_y;

			(*strain_xx_glo) = strain_xx;
			(*strain_yy_glo) = strain_yy;
			(*strain_xy_glo) = strain_xy;

			return i;
		}

		InverseMatrix_2D(temp_matrix);

		temp_dxi = temp_matrix[0][0] * (px - temp_x)
		           + temp_matrix[0][1] * (py - temp_y);
		temp_deta = temp_matrix[1][0] * (px - temp_x)
		            + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_deta;
		if (temp_xi < knot_vec_xi[0])
			temp_xi = knot_vec_xi[0];
		if (temp_xi > knot_vec_xi[cntl_p_n_xi + order_xi])
			temp_xi = knot_vec_xi[cntl_p_n_xi + order_xi];
		if (temp_eta < knot_vec_eta[0])
			temp_eta = knot_vec_eta[0];
		if (temp_eta > knot_vec_eta[cntl_p_n_eta + order_eta])
			temp_eta = knot_vec_eta[cntl_p_n_eta + order_eta];

		//temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
		//printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
	}
	temp_xi = knot_vec_xi[0] + knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = knot_vec_eta[0] + knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;
	for (i = 0; i < 1000; i++) {
		rlNURBS_surface(knot_vec_xi, knot_vec_eta,
		               cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
		               weight, order_xi, order_eta,
		               temp_xi, temp_eta,
		               &temp_x, &temp_y,
		               &temp_matrix[0][0], &temp_matrix[0][1],
		               &temp_matrix[1][0], &temp_matrix[1][1]);

		temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;
		
		//収束した場合////////////////////////////////////////////////////////////////
		//if (temp_tol_x < tol && temp_tol_y < tol) {
        if (temp_tol_x + temp_tol_y < tol) {
			printf("rlNURBS\n");
			if (temp_xi == knot_vec_xi[0])
			{
				break;
			}
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;

			int knot_n_xi = cntl_p_n_xi + order_xi + 1;
			int knot_n_eta = cntl_p_n_eta + order_eta + 1;
			double dtilda_xi = 0.0;
			double dtilda_eta = 0.0;
			double disp_x, disp_y;
			double dxi_x, deta_x, dxi_y, deta_y;
			double dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y;
			double temp_matrix2[2][2];
			double temp;
			double strain_xx, strain_yy, strain_xy;
			double stress_xx, stress_yy, stress_xy;

			for (i = 0; i < knot_n_xi; i++) {
				if ( knot_vec_xi[i] < temp_xi && temp_xi <= knot_vec_xi[i + 1]) {
					dtilda_xi = ( knot_vec_xi[i + 1] - knot_vec_xi[i] ) / 2.0;
					//printf("%f\n", dtilda_xi);
					break;
				}
			}
			for (i = 0; i < knot_n_eta; i++) {
				if ( knot_vec_eta[i] <= temp_eta && temp_eta < knot_vec_eta[i + 1]) {
					dtilda_eta = ( knot_vec_eta[i + 1] - knot_vec_eta[i] ) / 2.0;
					//printf("%f\n", dtilda_eta);
					break;
				}
			}

			rlNURBS_surface(knot_vec_xi, knot_vec_eta,
			               cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
			               weight, order_xi, order_eta,
			               temp_xi, temp_eta,
			               &temp, &temp,
			               &dxi_x, &deta_x,
			               &dxi_y, &deta_y);
			//printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
			//	   dxi_x, deta_x, dxi_y, deta_y);

			rlNURBS_surface(knot_vec_xi, knot_vec_eta,
			               disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
			               weight, order_xi, order_eta,
			               temp_xi, temp_eta,
			               &disp_x, &disp_y,
			               &dxi_disp_x, &deta_disp_x,
			               &dxi_disp_y, &deta_disp_y);
			//printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
			//	   dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y);

			temp_matrix2[0][0] = dxi_x * dtilda_xi;
			temp_matrix2[0][1] = dxi_y * dtilda_xi;
			temp_matrix2[1][0] = deta_x * dtilda_eta;
			temp_matrix2[1][1] = deta_y * dtilda_eta;

			InverseMatrix_2D(temp_matrix2);

			strain_xx = temp_matrix2[0][0] * dtilda_xi
			            * dxi_disp_x
			            + temp_matrix2[0][1] * dtilda_eta
			            * deta_disp_x;
			strain_yy = temp_matrix2[1][0] * dtilda_xi
			            * dxi_disp_y
			            + temp_matrix2[1][1] * dtilda_eta
			            * deta_disp_y;
			strain_xy = temp_matrix2[1][0] * dtilda_xi
			            * dxi_disp_x
			            + temp_matrix2[1][1] * dtilda_eta
			            * deta_disp_x
			            + temp_matrix2[0][0] * dtilda_xi
			            * dxi_disp_y
			            + temp_matrix2[0][1] * dtilda_eta
			            * deta_disp_y;

			double D_matrix[3][3] = {{0.0}};
			int stress_type_flag = 1;
			if (stress_type_flag == 0) { //平面応力状態
				temp = E * (1.0 - nu * nu);
				D_matrix[0][0] = temp;
				D_matrix[0][1] = nu * temp;
				D_matrix[1][0] = nu * temp;
				D_matrix[1][1] = temp;
				D_matrix[2][2] = (1.0 - nu) / 2.0 * temp;
			} else if (stress_type_flag == 1) { //平面ひずみ状態(2Dの場合はこっち)
				temp = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
				D_matrix[0][0] = temp;
				D_matrix[0][1] = nu / (1.0 - nu) * temp;
				D_matrix[1][0] = nu / (1.0 - nu) * temp;
				D_matrix[1][1] = temp;
				D_matrix[2][2] = (1.0 - 2.0 * nu) / 2.0 / (1.0 - nu) * temp;
			}

			stress_xx = D_matrix[0][0] * strain_xx
			            + D_matrix[0][1] * strain_yy;
			stress_yy = D_matrix[1][0] * strain_xx
			            + D_matrix[1][1] * strain_yy;
			stress_xy = D_matrix[2][2] * strain_xy;

			temp = sqrt(disp_x * disp_x + disp_y * disp_y);
			printf("x:   % 1.8e\n", px);
			printf("y:   % 1.8e\n", py);
			printf("xi:  % 1.8e\n", temp_xi);
			printf("eta: % 1.8e\n", temp_eta);
			printf("Displacement x: % 1.8e\n", disp_x);
			printf("Displacement y: % 1.8e\n", disp_y);
			printf("Displacement  : % 1.8e\n", temp);
			printf("Strain xx: % 1.8e\n", strain_xx);
			printf("Strain yy: % 1.8e\n", strain_yy);
			printf("Strain xy: % 1.8e\n", strain_xy);
			printf("Stress xx: % 1.8e\n", stress_xx);
			printf("Stress yy: % 1.8e\n", stress_yy);
			printf("Stress xy: % 1.8e\n", stress_xy);

			(*disp_x_glo) = disp_x;
			(*disp_y_glo) = disp_y;

			(*strain_xx_glo) = strain_xx;
			(*strain_yy_glo) = strain_yy;
			(*strain_xy_glo) = strain_xy;

			return i;
		}

		InverseMatrix_2D(temp_matrix);

		temp_dxi = temp_matrix[0][0] * (px - temp_x)
		           + temp_matrix[0][1] * (py - temp_y);
		temp_deta = temp_matrix[1][0] * (px - temp_x)
		            + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_deta;
		if (temp_xi < knot_vec_xi[0])
			temp_xi = knot_vec_xi[0];
		if (temp_xi > knot_vec_xi[cntl_p_n_xi + order_xi])
			temp_xi = knot_vec_xi[cntl_p_n_xi + order_xi];
		if (temp_eta < knot_vec_eta[0])
			temp_eta = knot_vec_eta[0];
		if (temp_eta > knot_vec_eta[cntl_p_n_eta + order_eta])
			temp_eta = knot_vec_eta[cntl_p_n_eta + order_eta];

		//temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
		//printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
	}
	temp_xi = knot_vec_xi[0] + knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = knot_vec_eta[0] + knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;
	for (i = 0; i < 1000; i++) {
		lrNURBS_surface(knot_vec_xi, knot_vec_eta,
		               cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
		               weight, order_xi, order_eta,
		               temp_xi, temp_eta,
		               &temp_x, &temp_y,
		               &temp_matrix[0][0], &temp_matrix[0][1],
		               &temp_matrix[1][0], &temp_matrix[1][1]);

		temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;
		
		//収束した場合////////////////////////////////////////////////////////////////
		//if (temp_tol_x < tol && temp_tol_y < tol) {
        if (temp_tol_x + temp_tol_y < tol) {
			printf("lrNURBS\n");
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;

			int knot_n_xi = cntl_p_n_xi + order_xi + 1;
			int knot_n_eta = cntl_p_n_eta + order_eta + 1;
			double dtilda_xi = 0.0;
			double dtilda_eta = 0.0;
			double disp_x, disp_y;
			double dxi_x, deta_x, dxi_y, deta_y;
			double dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y;
			double temp_matrix2[2][2];
			double temp;
			double strain_xx, strain_yy, strain_xy;
			double stress_xx, stress_yy, stress_xy;

			for (i = 0; i < knot_n_xi; i++) {
				if ( knot_vec_xi[i] <= temp_xi && temp_xi < knot_vec_xi[i + 1]) {
					dtilda_xi = ( knot_vec_xi[i + 1] - knot_vec_xi[i] ) / 2.0;
					//printf("%f\n", dtilda_xi);
					break;
				}
			}
			for (i = 0; i < knot_n_eta; i++) {
				if ( knot_vec_eta[i] < temp_eta && temp_eta <= knot_vec_eta[i + 1]) {
					dtilda_eta = ( knot_vec_eta[i + 1] - knot_vec_eta[i] ) / 2.0;
					//printf("%f\n", dtilda_eta);
					break;
				}
			}

			lrNURBS_surface(knot_vec_xi, knot_vec_eta,
			               cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
			               weight, order_xi, order_eta,
			               temp_xi, temp_eta,
			               &temp, &temp,
			               &dxi_x, &deta_x,
			               &dxi_y, &deta_y);
			//printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
			//	   dxi_x, deta_x, dxi_y, deta_y);

			lrNURBS_surface(knot_vec_xi, knot_vec_eta,
			               disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
			               weight, order_xi, order_eta,
			               temp_xi, temp_eta,
			               &disp_x, &disp_y,
			               &dxi_disp_x, &deta_disp_x,
			               &dxi_disp_y, &deta_disp_y);
			//printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
			//	   dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y);

			temp_matrix2[0][0] = dxi_x * dtilda_xi;
			temp_matrix2[0][1] = dxi_y * dtilda_xi;
			temp_matrix2[1][0] = deta_x * dtilda_eta;
			temp_matrix2[1][1] = deta_y * dtilda_eta;

			InverseMatrix_2D(temp_matrix2);

			strain_xx = temp_matrix2[0][0] * dtilda_xi
			            * dxi_disp_x
			            + temp_matrix2[0][1] * dtilda_eta
			            * deta_disp_x;
			strain_yy = temp_matrix2[1][0] * dtilda_xi
			            * dxi_disp_y
			            + temp_matrix2[1][1] * dtilda_eta
			            * deta_disp_y;
			strain_xy = temp_matrix2[1][0] * dtilda_xi
			            * dxi_disp_x
			            + temp_matrix2[1][1] * dtilda_eta
			            * deta_disp_x
			            + temp_matrix2[0][0] * dtilda_xi
			            * dxi_disp_y
			            + temp_matrix2[0][1] * dtilda_eta
			            * deta_disp_y;

			double D_matrix[3][3] = {{0.0}};
			int stress_type_flag = 1;
			if (stress_type_flag == 0) { //平面応力状態
				temp = E * (1.0 - nu * nu);
				D_matrix[0][0] = temp;
				D_matrix[0][1] = nu * temp;
				D_matrix[1][0] = nu * temp;
				D_matrix[1][1] = temp;
				D_matrix[2][2] = (1.0 - nu) / 2.0 * temp;
			} else if (stress_type_flag == 1) { //平面ひずみ状態(2Dの場合はこっち)
				temp = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
				D_matrix[0][0] = temp;
				D_matrix[0][1] = nu / (1.0 - nu) * temp;
				D_matrix[1][0] = nu / (1.0 - nu) * temp;
				D_matrix[1][1] = temp;
				D_matrix[2][2] = (1.0 - 2.0 * nu) / 2.0 / (1.0 - nu) * temp;
			}

			stress_xx = D_matrix[0][0] * strain_xx
			            + D_matrix[0][1] * strain_yy;
			stress_yy = D_matrix[1][0] * strain_xx
			            + D_matrix[1][1] * strain_yy;
			stress_xy = D_matrix[2][2] * strain_xy;

			temp = sqrt(disp_x * disp_x + disp_y * disp_y);
			printf("x:   % 1.8e\n", px);
			printf("y:   % 1.8e\n", py);
			printf("xi:  % 1.8e\n", temp_xi);
			printf("eta: % 1.8e\n", temp_eta);
			printf("Displacement x: % 1.8e\n", disp_x);
			printf("Displacement y: % 1.8e\n", disp_y);
			printf("Displacement  : % 1.8e\n", temp);
			printf("Strain xx: % 1.8e\n", strain_xx);
			printf("Strain yy: % 1.8e\n", strain_yy);
			printf("Strain xy: % 1.8e\n", strain_xy);
			printf("Stress xx: % 1.8e\n", stress_xx);
			printf("Stress yy: % 1.8e\n", stress_yy);
			printf("Stress xy: % 1.8e\n", stress_xy);

			(*disp_x_glo) = disp_x;
			(*disp_y_glo) = disp_y;

			(*strain_xx_glo) = strain_xx;
			(*strain_yy_glo) = strain_yy;
			(*strain_xy_glo) = strain_xy;

			return i;
		}

		InverseMatrix_2D(temp_matrix);

		temp_dxi = temp_matrix[0][0] * (px - temp_x)
		           + temp_matrix[0][1] * (py - temp_y);
		temp_deta = temp_matrix[1][0] * (px - temp_x)
		            + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_deta;
		if (temp_xi < knot_vec_xi[0])
			temp_xi = knot_vec_xi[0];
		if (temp_xi > knot_vec_xi[cntl_p_n_xi + order_xi])
			temp_xi = knot_vec_xi[cntl_p_n_xi + order_xi];
		if (temp_eta < knot_vec_eta[0])
			temp_eta = knot_vec_eta[0];
		if (temp_eta > knot_vec_eta[cntl_p_n_eta + order_eta])
			temp_eta = knot_vec_eta[cntl_p_n_eta + order_eta];

		//temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
		//printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
	}
	//printf("% 1.8e % 1.8e\n", temp_x, temp_y);
	return 0;
}

static void Calculation (char *filename_output,
						 int order_xi, int order_eta,
						 int knot_n_xi, int knot_n_eta,
						 int cntl_p_n_xi, int cntl_p_n_eta,
						 double *knot_vec_xi, double *knot_vec_eta,
						 double *cntl_px, double *cntl_py,
						 double *disp_cntl_px, double *disp_cntl_py,
						 double *weight) {
	int i, j, k, l;
	double temp1, temp2, temp3;
	double temp_matrix[2][2];

	//計算するξ,ηの値決定と ∂ξ/∂チルダξ, ∂η/∂チルダη の計算
	double calc_xi[MAX_POINTS];		//計算するξの値
	double calc_eta[MAX_POINTS];		//計算するηの値
	double dtilda_xi[MAX_KNOTS];		// ∂ξ/∂チルダξ
	double dtilda_eta[MAX_KNOTS];	// ∂η/∂チルダη
	k = 0;
	l = 0;
	for (i = 0; i < knot_n_xi - 1; i++) {
		if ( knot_vec_xi[i] != knot_vec_xi[i + 1] ) {
			calc_xi[k] = knot_vec_xi[i];
			printf("%d\t%f\n", k, calc_xi[k]);
			dtilda_xi[l] = ( knot_vec_xi[i + 1] - knot_vec_xi[i] ) / 2.0;
			printf("%d\t%f\n", k, dtilda_xi[k]);
			k++;
			l++;
			if (division_ele_xi > 1) {
				temp1 = (knot_vec_xi[i + 1] - knot_vec_xi[i])
						/ (double)division_ele_xi;
				for (j = 1; j < division_ele_xi; j++) {
					calc_xi[k] = calc_xi[k - 1] + temp1;
					printf("%d\t%f\n", k, calc_xi[k]);
					k++;
				}
			}
		}
	}
	calc_xi[k] = knot_vec_xi[knot_n_xi - 1];
	printf("%d\t%f\n", k, calc_xi[k]);
	//printf("\n");
	division_n_xi = k + 1;
	element_n_xi = l;

	k = 0;
	l = 0;
	for (i = 0; i < knot_n_eta - 1; i++) {
		if ( knot_vec_eta[i] != knot_vec_eta[i + 1] ) {
			calc_eta[k] = knot_vec_eta[i];
			//printf("%d\t%f\n", k, calc_eta[k]);
			dtilda_eta[l] = ( knot_vec_eta[i + 1] - knot_vec_eta[i] ) / 2.0;
			//printf("%d\t%f\n", k, dtilda_eta[k]);
			k++;
			l++;
			if (division_ele_eta > 1) {
				temp1 = (knot_vec_eta[i + 1] - knot_vec_eta[i])
						/ (double)division_ele_eta;
				for (j = 1; j < division_ele_eta; j++) {
					calc_eta[k] = calc_eta[k - 1] + temp1;
					//printf("%d\t%f\n", k, calc_eta[k]);
					k++;
				}
			}
		}
	}
	calc_eta[k] = knot_vec_eta[knot_n_eta - 1];
	//printf("%d\t%f\n", k, calc_eta[k]);
	//printf("\n");
	division_n_eta = k + 1;
	element_n_eta = l;

	if (element_n_xi > MAX_ELEMENTS) {
		printf("Error!!\n");
		printf("Too many elements at xi!\n"
			   "Maximum of elements is %d (Now %d)\n"
			   "\n", MAX_ELEMENTS, element_n_xi);
		exit(1);
	}
	if (element_n_eta > MAX_ELEMENTS) {
		printf("Error!!\n");
		printf("Too many elements at eta!\n"
			   "Maximum of elements is %d (Now %d)\n"
			   "\n", MAX_ELEMENTS, element_n_eta);
		exit(1);
	}

	int ii, jj, kk, ll;

	//メッシュ座標計算
	printf("Start Calculation mesh\n\n");
	for (i = 0; i < division_n_xi; i++) {
		ii = i / division_ele_xi;
		kk = i % division_ele_xi;
		for (j = 0; j < division_n_eta; j++) {
			jj = j / division_ele_eta;
			ll = j % division_ele_eta;
			lNURBS_surface(knot_vec_xi, knot_vec_eta,
						   cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
						   weight, order_xi, order_eta,
						   calc_xi[i], calc_eta[j],
						   &coord_x[i][j], &coord_y[i][j],
						   &dxi_x[ii][jj][kk][ll], &deta_x[ii][jj][kk][ll],
						   &dxi_y[ii][jj][kk][ll], &deta_y[ii][jj][kk][ll]);
			printf("[%d][%d] [%d][%d][%d][%d]"
				   "% 1.4e % 1.4e "
				   "% 1.4e % 1.4e\n",
				   i, j, ii, jj, kk, ll,
				   calc_xi[i], calc_eta[j],
				   coord_x[i][j], coord_y[i][j]);
		}
		//printf("\n");
	}
	printf("\n");
	printf("End Calculation mesh\n\n");

	if (fields_flag) {
		//変位計算
		printf("Start Calculation displpacement\n\n");
		for (i = 0; i < division_n_xi; i++) {
			ii = i / division_ele_xi;
			kk = i % division_ele_xi;
			for (j = 0; j < division_n_eta; j++) {
				jj = j / division_ele_eta;
				ll = j % division_ele_eta;
				lNURBS_surface(knot_vec_xi, knot_vec_eta,
							   disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
							   weight, order_xi, order_eta,
							   calc_xi[i], calc_eta[j],
							   &disp_x[i][j], &disp_y[i][j],
							   &dxi_disp_x[ii][jj][kk][ll], &deta_disp_x[ii][jj][kk][ll],
							   &dxi_disp_y[ii][jj][kk][ll], &deta_disp_y[ii][jj][kk][ll]);
				printf("[%d][%d] [%d][%d][%d][%d]"
					   "% 1.4e % 1.4e "
					   "% 1.4e % 1.4e\n",
					   i, j, ii, jj, kk, ll,
					   calc_xi[i], calc_eta[j],
					   disp_x[i][j], disp_y[i][j]);
			}
			//printf("\n");
		}
		printf("\n");
		printf("End Calculation displpacement\n\n");

		//足りない微分値計算
		for (ii = 0; ii < element_n_xi; ii++) {
			for (jj = 0; jj < element_n_eta; jj++) {
				kk = division_ele_xi;
				i = (ii + 1) * division_ele_xi;
				j = jj * division_ele_eta;
				for (ll = 1; ll < division_ele_eta; ll++) {
					j++;
					rNURBS_surface(knot_vec_xi, knot_vec_eta,
								   cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
								   weight, order_xi, order_eta,
								   calc_xi[i], calc_eta[j],
								   &coord_x[i][j], &coord_y[i][j],
								   &dxi_x[ii][jj][kk][ll], &deta_x[ii][jj][kk][ll],
								   &dxi_y[ii][jj][kk][ll], &deta_y[ii][jj][kk][ll]);
					rNURBS_surface(knot_vec_xi, knot_vec_eta,
								   disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
								   weight, order_xi, order_eta,
								   calc_xi[i], calc_eta[j],
								   &disp_x[i][j], &disp_y[i][j],
								   &dxi_disp_x[ii][jj][kk][ll], &deta_disp_x[ii][jj][kk][ll],
								   &dxi_disp_y[ii][jj][kk][ll], &deta_disp_y[ii][jj][kk][ll]);
					/*
					printf("[%d][%d] [%d][%d][%d][%d]"
						   "% 1.4e % 1.4e % 1.4e % 1.4e"
						   "% 1.4e % 1.4e % 1.4e % 1.4e\n",
						   i, j, ii, jj, kk, ll,
						   dxi_x[ii][jj][kk][ll], deta_x[ii][jj][kk][ll],
						   dxi_y[ii][jj][kk][ll], deta_y[ii][jj][kk][ll],
						   dxi_disp_x[ii][jj][kk][ll], deta_disp_x[ii][jj][kk][ll],
						   dxi_disp_y[ii][jj][kk][ll], deta_disp_y[ii][jj][kk][ll]);
					*/
				}

				ll = division_ele_eta;
				i = ii * division_ele_xi;
				j = (jj + 1) * division_ele_eta;
				for (kk = 1; kk <= division_ele_xi; kk++) {
					i++;
					rNURBS_surface(knot_vec_xi, knot_vec_eta,
								   cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
								   weight, order_xi, order_eta,
								   calc_xi[i], calc_eta[j],
								   &coord_x[i][j], &coord_y[i][j],
								   &dxi_x[ii][jj][kk][ll], &deta_x[ii][jj][kk][ll],
								   &dxi_y[ii][jj][kk][ll], &deta_y[ii][jj][kk][ll]);
					rNURBS_surface(knot_vec_xi, knot_vec_eta,
								   disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
								   weight, order_xi, order_eta,
								   calc_xi[i], calc_eta[j],
								   &disp_x[i][j], &disp_y[i][j],
								   &dxi_disp_x[ii][jj][kk][ll], &deta_disp_x[ii][jj][kk][ll],
								   &dxi_disp_y[ii][jj][kk][ll], &deta_disp_y[ii][jj][kk][ll]);
					/*
					printf("[%d][%d] [%d][%d][%d][%d]"
						   "% 1.4e % 1.4e % 1.4e % 1.4e"
						   "% 1.4e % 1.4e % 1.4e % 1.4e\n",
						   i, j, ii, jj, kk, ll,
						   dxi_x[ii][jj][kk][ll], deta_x[ii][jj][kk][ll],
						   dxi_y[ii][jj][kk][ll], deta_y[ii][jj][kk][ll],
						   dxi_disp_x[ii][jj][kk][ll], deta_disp_x[ii][jj][kk][ll],
						   dxi_disp_y[ii][jj][kk][ll], deta_disp_y[ii][jj][kk][ll]);
					*/
				}

				kk = division_ele_xi;
				ll = 0;
				i = (ii + 1) * division_ele_xi;
				j = jj * division_ele_eta;
				rlNURBS_surface(knot_vec_xi, knot_vec_eta,
								cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
								weight, order_xi, order_eta,
								calc_xi[i], calc_eta[j],
								&coord_x[i][j], &coord_y[i][j],
								&dxi_x[ii][jj][kk][ll], &deta_x[ii][jj][kk][ll],
								&dxi_y[ii][jj][kk][ll], &deta_y[ii][jj][kk][ll]);
				rlNURBS_surface(knot_vec_xi, knot_vec_eta,
								disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
								weight, order_xi, order_eta,
								calc_xi[i], calc_eta[j],
								&disp_x[i][j], &disp_y[i][j],
								&dxi_disp_x[ii][jj][kk][ll], &deta_disp_x[ii][jj][kk][ll],
								&dxi_disp_y[ii][jj][kk][ll], &deta_disp_y[ii][jj][kk][ll]);
				/*
					printf("[%d][%d] [%d][%d][%d][%d]"
						   "% 1.4e % 1.4e % 1.4e % 1.4e"
						   "% 1.4e % 1.4e % 1.4e % 1.4e\n",
						   i, j, ii, jj, kk, ll,
						   dxi_x[ii][jj][kk][ll], deta_x[ii][jj][kk][ll],
						   dxi_y[ii][jj][kk][ll], deta_y[ii][jj][kk][ll],
						   dxi_disp_x[ii][jj][kk][ll], deta_disp_x[ii][jj][kk][ll],
						   dxi_disp_y[ii][jj][kk][ll], deta_disp_y[ii][jj][kk][ll]);
				*/

				kk = 0;
				ll = division_ele_eta;
				i = ii * division_ele_xi;
				j = (jj + 1) * division_ele_eta;
				lrNURBS_surface(knot_vec_xi, knot_vec_eta,
								cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
								weight, order_xi, order_eta,
								calc_xi[i], calc_eta[j],
								&coord_x[i][j], &coord_y[i][j],
								&dxi_x[ii][jj][kk][ll], &deta_x[ii][jj][kk][ll],
								&dxi_y[ii][jj][kk][ll], &deta_y[ii][jj][kk][ll]);
				lrNURBS_surface(knot_vec_xi, knot_vec_eta,
								disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
								weight, order_xi, order_eta,
								calc_xi[i], calc_eta[j],
								&disp_x[i][j], &disp_y[i][j],
								&dxi_disp_x[ii][jj][kk][ll], &deta_disp_x[ii][jj][kk][ll],
								&dxi_disp_y[ii][jj][kk][ll], &deta_disp_y[ii][jj][kk][ll]);
				/*
					printf("[%d][%d] [%d][%d][%d][%d]"
						   "% 1.4e % 1.4e % 1.4e % 1.4e"
						   "% 1.4e % 1.4e % 1.4e % 1.4e\n",
						   i, j, ii, jj, kk, ll,
						   dxi_x[ii][jj][kk][ll], deta_x[ii][jj][kk][ll],
						   dxi_y[ii][jj][kk][ll], deta_y[ii][jj][kk][ll],
						   dxi_disp_x[ii][jj][kk][ll], deta_disp_x[ii][jj][kk][ll],
						   dxi_disp_y[ii][jj][kk][ll], deta_disp_y[ii][jj][kk][ll]);
				printf("\n");
				*/
			}
		}

		/*
		for (ii = 0; ii < element_n_xi; ii++) {
			for (jj = 0; jj < element_n_eta; jj++) {
				for (kk = 0; kk <= division_ele_xi; kk++) {
					for (ll = 0; ll <= division_ele_eta; ll++) {
						printf("[%d][%d][%d][%d]"
							   "% 1.4e % 1.4e % 1.4e % 1.4e"
							   "% 1.4e % 1.4e % 1.4e % 1.4e\n",
							   ii, jj, kk, ll,
							   dxi_x[ii][jj][kk][ll], deta_x[ii][jj][kk][ll],
							   dxi_y[ii][jj][kk][ll], deta_y[ii][jj][kk][ll],
							   dxi_disp_x[ii][jj][kk][ll], deta_disp_x[ii][jj][kk][ll],
							   dxi_disp_y[ii][jj][kk][ll], deta_disp_y[ii][jj][kk][ll]);
					}
				}
			}
		}
		printf("\n");
		*/

		//ひずみ計算
		printf("Start Calculation Strain\n\n");
		for (i = 0; i < element_n_xi; i++) {
			for (j = 0; j < element_n_eta; j++) {
				temp1 = dtilda_xi[i];
				temp2 = dtilda_eta[j];
				for (k = 0; k < division_ele_xi + 1; k++) {
					for (l = 0; l < division_ele_eta + 1; l++) {
						temp_matrix[0][0] = dxi_x[i][j][k][l] * temp1;
						temp_matrix[0][1] = dxi_y[i][j][k][l] * temp1;
						temp_matrix[1][0] = deta_x[i][j][k][l] * temp2;
						temp_matrix[1][1] = deta_y[i][j][k][l] * temp2;

						InverseMatrix_2D(temp_matrix);

						strain_xx[i][j][k][l] = temp_matrix[0][0] * temp1
												* dxi_disp_x[i][j][k][l]
												+ temp_matrix[0][1] * temp2
												* deta_disp_x[i][j][k][l];
						strain_yy[i][j][k][l] = temp_matrix[1][0] * temp1
												* dxi_disp_y[i][j][k][l]
												+ temp_matrix[1][1] * temp2
												* deta_disp_y[i][j][k][l];
						strain_xy[i][j][k][l] = temp_matrix[1][0] * temp1
												* dxi_disp_x[i][j][k][l]
												+ temp_matrix[1][1] * temp2
												* deta_disp_x[i][j][k][l]
												+ temp_matrix[0][0] * temp1
												* dxi_disp_y[i][j][k][l]
												+ temp_matrix[0][1] * temp2
												* deta_disp_y[i][j][k][l];

						printf("[%d][%d][%d][%d]\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n",
							   i, j, k, l, temp1, temp2,
							   strain_xx[i][j][k][l], strain_yy[i][j][k][l], strain_xy[i][j][k][l]);
					}
				}
			}
			printf("\n");
		}
		printf("End Calculation Strain\n\n");

		//Dマトリクスの計算
		double D_matrix[3][3] = {{0.0}};
		if (stress_type_flag == 0) { //平面応力状態
			temp1 = E * (1.0 - nu * nu);
			D_matrix[0][0] = temp1;
			D_matrix[0][1] = nu * temp1;
			D_matrix[1][0] = nu * temp1;
			D_matrix[1][1] = temp1;
			D_matrix[2][2] = (1.0 - nu) / 2.0 * temp1;
		} else if (stress_type_flag == 1) { //平面ひずみ状態(2Dの場合はこっち)
			temp1 = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
			D_matrix[0][0] = temp1;
			D_matrix[0][1] = nu / (1.0 - nu) * temp1;
			D_matrix[1][0] = nu / (1.0 - nu) * temp1;
			D_matrix[1][1] = temp1;
			D_matrix[2][2] = (1.0 - 2.0 * nu) / 2.0 / (1.0 - nu) * temp1;
		}

		printf("Start Calculation Stress\n\n");
		for (i = 0; i < element_n_xi; i++) {
			for (j = 0; j < element_n_eta; j++) {
				for (k = 0; k < division_ele_xi + 1; k++) {
					for (l = 0; l < division_ele_eta + 1; l++) {
						stress_xx[i][j][k][l] = D_matrix[0][0] * strain_xx[i][j][k][l]
												+ D_matrix[0][1] * strain_yy[i][j][k][l];
						stress_yy[i][j][k][l] = D_matrix[1][0] * strain_xx[i][j][k][l]
												+ D_matrix[1][1] * strain_yy[i][j][k][l];
						stress_xy[i][j][k][l] = D_matrix[2][2] * strain_xy[i][j][k][l];
						printf("[%d][%d][%d][%d]\t% 1.4e\t% 1.4e\t% 1.4e\n", i, j, k, l,
							   stress_xx[i][j][k][l], stress_yy[i][j][k][l], stress_xy[i][j][k][l]);
					}
				}
				printf("\n");
			}
		}
		printf("End Calculation Stress\n\n");
	}

	//書き込み
	fp = fopen(filename_output, "a");
	if (fields_flag) {
		fprintf(fp, "%d\t%d\t%d\t%d\n",
				division_n_xi, division_n_eta,
				element_n_xi, element_n_eta);
		for (i = 0; i < division_n_xi; i++) {
			for (j = 0; j < division_n_eta; j++) {
				temp1 = disp_x[i][j];
				temp2 = disp_y[i][j];
				temp3 = sqrt(temp1 * temp1 + temp2 * temp2);
				fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n",
						coord_x[i][j], coord_y[i][j],
						disp_x[i][j], disp_y[i][j], temp3);
			}
		}
		for (i = 0; i < element_n_xi; i++) {
			for (j = 0; j < element_n_eta; j++) {
				for (k = 0; k < division_ele_xi + 1; k++) {
					for (l = 0; l < division_ele_eta + 1; l++) {
						fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n",
								strain_xx[i][j][k][l], strain_yy[i][j][k][l], strain_xy[i][j][k][l],
								stress_xx[i][j][k][l], stress_yy[i][j][k][l], stress_xy[i][j][k][l]);
					}
				}
			}
		}
	} else {
		fprintf(fp, "%d\t%d\t%d\t%d\n",
				division_n_xi, division_n_eta,
				element_n_xi, element_n_eta);
		for (i = 0; i < division_n_xi; i++) {
			for (j = 0; j < division_n_eta; j++) {
				fprintf(fp, "% 1.4e\t% 1.4e\n",
						coord_x[i][j], coord_y[i][j]);
			}
		}
	}
	fclose(fp);

	//グラフ用ファイル書き込み
	fp = fopen("disp_graph.txt", "a");
	for (i = 0; i < division_n_xi; i++) {
		for (j = 0; j < division_n_eta; j++) {
			temp1 = disp_x[i][j];
			temp2 = disp_y[i][j];
			//temp3 = sqrt(temp1 * temp1 + temp2 * temp2);
			fprintf(fp, "%d\t% 1.15e\t% 1.15e\t% 1.15e\t% 1.15e\n",
					graph_patch_n,
					coord_x[i][j], coord_y[i][j],
					disp_x[i][j], disp_y[i][j]);
		}
	}
	fclose(fp);
	
	fp = fopen("stress_y_graph.txt", "a");
 	for (i = 0; i < element_n_xi; i++) {
		for (j = 0; j < element_n_eta; j++) {
			for (k = 0; k < division_ele_xi + 1; k++) {
				for (l = 0; l < division_ele_eta + 1; l++) {
				fprintf(fp, "%d\t% 1.15e\t% 1.15e\t% 1.15e\n", 
							graph_patch_n,
							coord_x[i * division_ele_xi + k][j * division_ele_eta + l], 
							coord_y[i * division_ele_xi + k][j * division_ele_eta + l], 
							stress_yy[i][j][k][l]);
				}
			}
		}
	}
	fclose(fp);

	fp = fopen("stress_vm_graph.txt", "a");
 	for (i = 0; i < element_n_xi; i++) {
		for (j = 0; j < element_n_eta; j++) {
			for (k = 0; k < division_ele_xi + 1; k++) {
				for (l = 0; l < division_ele_eta + 1; l++) {
					double stress_vm;
					double sum = stress_xx[i][j][k][l] + stress_yy[i][j][k][l];
					double dif = stress_xx[i][j][k][l] - stress_yy[i][j][k][l];
					double tau2 = stress_xy[i][j][k][l] * stress_xy[i][j][k][l];
					double temp1, temp2;
					double temp3;
					temp1 = 0.5 * sum;
					temp2 = 0.5 * sqrt(dif * dif + 4 * tau2);
					stress_vm = sqrt(temp1 * temp1 + 3 * temp2 * temp2);
					temp3 = sqrt(calc_xi[i * division_ele_xi + k] * calc_xi[i * division_ele_xi + k]
								+ calc_eta[j * division_ele_eta + l] * calc_eta[j * division_ele_eta + l]);
					fprintf(fp, "% 1.10e\t% 1.10e\t% 1.10e\t% 1.10e\t"
								"% 1.10e\t% 1.10e\n", 
							calc_xi[i * division_ele_xi + k], 
							calc_eta[j * division_ele_eta + l],
							coord_x[i * division_ele_xi + k][j * division_ele_eta + l], 
							coord_y[i * division_ele_xi + k][j * division_ele_eta + l], 
							stress_vm, temp3);
				}
			}
		}
	}
	fclose(fp);
}

//for s-IGA
//重ね合わせた結果の出力
static void Calculation_overlay(char *filename_output_over,
						        int order_xi_loc, int order_eta_loc,
						        int knot_n_xi_loc, int knot_n_eta_loc,
					    	    int cntl_p_n_xi_loc, int cntl_p_n_eta_loc,
						        double *knot_vec_xi_loc, double *knot_vec_eta_loc,
						        double *cntl_px_loc, double *cntl_py_loc,
						        double *disp_cntl_px_loc, double *disp_cntl_py_loc,
						        double *weight_loc,
                                int order_xi_glo, int order_eta_glo,
						        int knot_n_xi_glo, int knot_n_eta_glo,
					    	    int cntl_p_n_xi_glo, int cntl_p_n_eta_glo,
						        double *knot_vec_xi_glo, double *knot_vec_eta_glo,
						        double *cntl_px_glo, double *cntl_py_glo,
						        double *disp_cntl_px_glo, double *disp_cntl_py_glo,
						        double *weight_glo)
{
	int i, j, k, l;
	double temp1, temp2, temp3;
	//double temp_matrix[2][2];

    double output_xi, output_eta;
	double disp_x_glo;
	double disp_y_glo;
    double strain_xx_glo = 0;
    double strain_yy_glo = 0;
    double strain_xy_glo = 0;
	//, strain_yy_glo, strain_xy_glo;

	//計算するξ,ηの値決定と ∂ξ/∂チルダξ, ∂η/∂チルダη の計算
	double calc_xi_loc[MAX_POINTS];		//計算するξの値local
	double calc_eta_loc[MAX_POINTS];		//計算するηの値local
	//double dtilda_xi[MAX_KNOTS];		// ∂ξ/∂チルダξ
	//double dtilda_eta[MAX_KNOTS];	// ∂η/∂チルダη
	k = 0;
	l = 0;
	for (i = 0; i < knot_n_xi_loc - 1; i++) {
		if ( knot_vec_xi_loc[i] != knot_vec_xi_loc[i + 1] ) {
			calc_xi_loc[k] = knot_vec_xi_loc[i];
			printf("%d\t%f\n", k, calc_xi_loc[k]);
			//dtilda_xi[l] = ( knot_vec_xi_loc[i + 1] - knot_vec_xi_loc[i] ) / 2.0;
			//printf("%d\t%f\n", k, dtilda_xi[k]);
			k++;
			l++;
			if (division_ele_xi > 1) {
				temp1 = (knot_vec_xi_loc[i + 1] - knot_vec_xi_loc[i])
						/ (double)division_ele_xi;
				for (j = 1; j < division_ele_xi; j++) {
					calc_xi_loc[k] = calc_xi_loc[k - 1] + temp1;
					printf("%d\t%f\n", k, calc_xi_loc[k]);
					k++;
				}
			}
		}
	}
	calc_xi_loc[k] = knot_vec_xi_loc[knot_n_xi_loc - 1];
	printf("%d\t%f\n", k, calc_xi_loc[k]);
	//printf("\n");
	division_n_xi = k + 1;
	element_n_xi = l;

	k = 0;
	l = 0;
	for (i = 0; i < knot_n_eta_loc - 1; i++) {
		if ( knot_vec_eta_loc[i] != knot_vec_eta_loc[i + 1] ) {
			calc_eta_loc[k] = knot_vec_eta_loc[i];
			//printf("%d\t%f\n", k, calc_eta[k]);
			//dtilda_eta[l] = ( knot_vec_eta_loc[i + 1] - knot_vec_eta_loc[i] ) / 2.0;
			//printf("%d\t%f\n", k, dtilda_eta[k]);
			k++;
			l++;
			if (division_ele_eta > 1) {
				temp1 = (knot_vec_eta_loc[i + 1] - knot_vec_eta_loc[i])
						/ (double)division_ele_eta;
				for (j = 1; j < division_ele_eta; j++) {
					calc_eta_loc[k] = calc_eta_loc[k - 1] + temp1;
					//printf("%d\t%f\n", k, calc_eta[k]);
					k++;
				}
			}
		}
	}
	calc_eta_loc[k] = knot_vec_eta_loc[knot_n_eta_loc - 1];
	//printf("%d\t%f\n", k, calc_eta[k]);
	//printf("\n");
	division_n_eta = k + 1;
	element_n_eta = l;

	if (element_n_xi > MAX_ELEMENTS) {
		printf("Error!!\n");
		printf("Too many elements at xi!\n"
			   "Maximum of elements is %d (Now %d)\n"
			   "\n", MAX_ELEMENTS, element_n_xi);
		exit(1);
	}
	if (element_n_eta > MAX_ELEMENTS) {
		printf("Error!!\n");
		printf("Too many elements at eta!\n"
			   "Maximum of elements is %d (Now %d)\n"
			   "\n", MAX_ELEMENTS, element_n_eta);
		exit(1);
	}

	int ii, jj, kk, ll;

	//メッシュ座標計算
	printf("Start Calculation overlay mesh\n\n");
	printf("Start Calculation overlay displpacement\n\n");
	printf("Start Calculation overlay Strain\n\n");
	for (i = 0; i < division_n_xi; i++) {
		ii = i / division_ele_xi;
		kk = i % division_ele_xi;
		for (j = 0; j < division_n_eta; j++) {
			jj = j / division_ele_eta;
			ll = j % division_ele_eta;
			lNURBS_surface(knot_vec_xi_loc, knot_vec_eta_loc,
						   cntl_px_loc, cntl_py_loc, cntl_p_n_xi_loc, cntl_p_n_eta_loc,
						   weight_loc, order_xi_loc, order_eta_loc,
						   calc_xi_loc[i], calc_eta_loc[j],
						   &coord_x[i][j], &coord_y[i][j],
						   &dxi_x[ii][jj][kk][ll], &deta_x[ii][jj][kk][ll],
						   &dxi_y[ii][jj][kk][ll], &deta_y[ii][jj][kk][ll]);
			printf("[%d][%d] [%d][%d][%d][%d]"
				   "% 1.4e % 1.4e "
				   "% 1.4e % 1.4e\n",
				   i, j, ii, jj, kk, ll,
				   calc_xi_loc[i], calc_eta_loc[j],
				   coord_x[i][j], coord_y[i][j]);
                        
            int itr_n = CalcXiEtaByNR(coord_x[i][j], coord_y[i][j],
                                      knot_vec_xi_glo, knot_vec_eta_glo,
				                      cntl_px_glo, cntl_py_glo,
			                          disp_cntl_px_glo, disp_cntl_py_glo,
			                          cntl_p_n_xi_glo, cntl_p_n_eta_glo,
			                          weight_glo, order_xi_glo, order_eta_glo,
			                          &output_xi, &output_eta,
									  &disp_x_glo, &disp_y_glo,
                                      &strain_xx_glo, &strain_yy_glo, &strain_xy_glo);
            if (itr_n == 0)
            {
                printf("itr=0\n");
            }
			printf("iteration : %d\n",itr_n);

			//ローカル内の表示点上のグローバル変位
			printf("disp_x_glo =% 1.4e\tdisp_y_glo =% 1.4e\n",
					disp_x_glo, disp_y_glo);
			printf("%1.4e\t%1.4e\n",disp_x[i][j],disp_y[i][j]);
			disp_x[i][j] += disp_x_glo;
			disp_y[i][j] += disp_y_glo;
			printf("% 1.4e\t% 1.4e\n",disp_x[i][j],disp_y[i][j]);

			//ローカル内の表示点上のグローバルひずみ
            printf("strain_xx_glo =% 1.4e\n"
                   "strain_yy_glo =% 1.4e\n"
                   "strain_xy_glo =% 1.4e\n",
                   strain_xx_glo, strain_yy_glo, strain_xy_glo);
			strain_xx[ii][jj][kk][ll] += strain_xx_glo;
			strain_yy[ii][jj][kk][ll] += strain_yy_glo;
			strain_xy[ii][jj][kk][ll] += strain_xy_glo;
			printf("test[%d][%d][%d][%d]\n",ii,jj,kk,ll);
			if (jj > 0 && ll == 0)
			{
				strain_xx[ii][jj-1][kk][division_ele_eta] += strain_xx_glo;
				strain_yy[ii][jj-1][kk][division_ele_eta] += strain_yy_glo;
				strain_xy[ii][jj-1][kk][division_ele_eta] += strain_xy_glo;
				printf("test[%d][%d][%d][%d]\n",ii,jj-1,kk,division_ele_eta);
			}
			if (ii > 0 && kk == 0)
			{
				strain_xx[ii-1][jj][division_ele_xi][ll] += strain_xx_glo;
				strain_yy[ii-1][jj][division_ele_xi][ll] += strain_yy_glo;
				strain_xy[ii-1][jj][division_ele_xi][ll] += strain_xy_glo;
				printf("test[%d][%d][%d][%d]\n",ii-1,jj,division_ele_xi,ll);
			}
			if (ii > 0 && jj > 0 && kk == 0 && ll == 0)
			{
				strain_xx[ii-1][jj-1][division_ele_xi][division_ele_eta] += strain_xx_glo;
				strain_yy[ii-1][jj-1][division_ele_xi][division_ele_eta] += strain_yy_glo;
				strain_xy[ii-1][jj-1][division_ele_xi][division_ele_eta] += strain_xy_glo;
				printf("test[%d][%d][%d][%d]\n",ii-1,jj-1,division_ele_xi,division_ele_eta);
			}
			printf("% 1.4e\t% 1.4e\t% 1.4e\n",
					strain_xx[ii][jj][kk][ll],
					strain_yy[ii][jj][kk][ll],
					strain_xy[ii][jj][kk][ll]);
        
		}
		//printf("\n");
	}
	printf("\n");
	printf("End Calculation overlay mesh\n\n");
	printf("End Calculation overlay displpacement\n\n");
	printf("End Calculation overlay Strain\n\n");
    
	//Dマトリクスの計算
	double D_matrix[3][3] = {{0.0}};
	if (stress_type_flag == 0) { //平面応力状態
		temp1 = E * (1.0 - nu * nu);
		D_matrix[0][0] = temp1;
		D_matrix[0][1] = nu * temp1;
		D_matrix[1][0] = nu * temp1;
		D_matrix[1][1] = temp1;
		D_matrix[2][2] = (1.0 - nu) / 2.0 * temp1;
	} else if (stress_type_flag == 1) { //平面ひずみ状態(2Dの場合はこっち)
		temp1 = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
		D_matrix[0][0] = temp1;
		D_matrix[0][1] = nu / (1.0 - nu) * temp1;
		D_matrix[1][0] = nu / (1.0 - nu) * temp1;
		D_matrix[1][1] = temp1;
		D_matrix[2][2] = (1.0 - 2.0 * nu) / 2.0 / (1.0 - nu) * temp1;
	}
	
	printf("Start Calculation overlay Stress\n\n");
	for (i = 0; i < element_n_xi; i++) {
		for (j = 0; j < element_n_eta; j++) {
			for (k = 0; k < division_ele_xi + 1; k++) {
				for (l = 0; l < division_ele_eta + 1; l++) {
					stress_xx[i][j][k][l] = D_matrix[0][0] * strain_xx[i][j][k][l]
											+ D_matrix[0][1] * strain_yy[i][j][k][l];
					stress_yy[i][j][k][l] = D_matrix[1][0] * strain_xx[i][j][k][l]
											+ D_matrix[1][1] * strain_yy[i][j][k][l];
					stress_xy[i][j][k][l] = D_matrix[2][2] * strain_xy[i][j][k][l];
					printf("[%d][%d][%d][%d]\t% 1.4e\t% 1.4e\t% 1.4e\n", i, j, k, l,
						   stress_xx[i][j][k][l], stress_yy[i][j][k][l], stress_xy[i][j][k][l]);
				}
			}
			printf("\n");
		}
	}
	printf("End Calculation overlay Stress\n\n");	



	//書き込み
	fp = fopen(filename_output_over, "a");
	if (fields_flag) {
		fprintf(fp, "%d\t%d\t%d\t%d\n",
				division_n_xi, division_n_eta,
				element_n_xi, element_n_eta);
		for (i = 0; i < division_n_xi; i++) {
			for (j = 0; j < division_n_eta; j++) {
				temp1 = disp_x[i][j];
				temp2 = disp_y[i][j];
				temp3 = sqrt(temp1 * temp1 + temp2 * temp2);
				fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n",
						coord_x[i][j], coord_y[i][j],
						disp_x[i][j], disp_y[i][j], temp3);
			}
		}
		for (i = 0; i < element_n_xi; i++) {
			for (j = 0; j < element_n_eta; j++) {
				for (k = 0; k < division_ele_xi + 1; k++) {
					for (l = 0; l < division_ele_eta + 1; l++) {
						fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n",
								strain_xx[i][j][k][l], strain_yy[i][j][k][l], strain_xy[i][j][k][l],
								stress_xx[i][j][k][l], stress_yy[i][j][k][l], stress_xy[i][j][k][l]);
					}
				}
			}
		}
	} else {
		fprintf(fp, "%d\t%d\t%d\t%d\n",
				division_n_xi, division_n_eta,
				element_n_xi, element_n_eta);
		for (i = 0; i < division_n_xi; i++) {
			for (j = 0; j < division_n_eta; j++) {
				fprintf(fp, "% 1.4e\t% 1.4e\n",
						coord_x[i][j], coord_y[i][j]);
			}
		}
	}
	fclose(fp);

	//グラフ用ファイル書き込み
	fp = fopen("over_disp_graph.txt", "a");
	for (i = 0; i < division_n_xi; i++) {
		for (j = 0; j < division_n_eta; j++) {
			temp1 = disp_x[i][j];
			temp2 = disp_y[i][j];
			//temp3 = sqrt(temp1 * temp1 + temp2 * temp2);
			fprintf(fp, "%d\t% 1.15e\t% 1.15e\t% 1.15e\t% 1.15e\n",
					graph_patch_n,
					coord_x[i][j], coord_y[i][j],
					disp_x[i][j], disp_y[i][j]);
		}
	}
	fclose(fp);

	fp = fopen("over_stress_x_graph.txt", "a");
 	for (i = 0; i < element_n_xi; i++) {
		for (j = 0; j < element_n_eta; j++) {
			for (k = 0; k < division_ele_xi + 1; k++) {
				for (l = 0; l < division_ele_eta + 1; l++) {
				fprintf(fp, "% 1.10e\t% 1.10e\t% 1.10e\n", 
							coord_x[i * division_ele_xi + k][j * division_ele_eta + l], 
							coord_y[i * division_ele_xi + k][j * division_ele_eta + l], 
							stress_xx[i][j][k][l]);
				}
			}
		}
	}
	fclose(fp);

	fp = fopen("over_stress_y_graph.txt", "a");
 	for (i = 0; i < element_n_xi; i++) {
		for (j = 0; j < element_n_eta; j++) {
			for (k = 0; k < division_ele_xi + 1; k++) {
				for (l = 0; l < division_ele_eta + 1; l++) {
				fprintf(fp, "% 1.15e\t% 1.15e\t% 1.15e\n", 
							coord_x[i * division_ele_xi + k][j * division_ele_eta + l], 
							coord_y[i * division_ele_xi + k][j * division_ele_eta + l], 
							stress_yy[i][j][k][l]);
				}
			}
		}
	}
	fclose(fp);

	fp = fopen("over_stress_r_graph.txt", "a");
 	for (i = 0; i < element_n_xi; i++) {
		for (j = 0; j < element_n_eta; j++) {
			for (k = 0; k < division_ele_xi + 1; k++) {
				for (l = 0; l < division_ele_eta + 1; l++) {
					double stress_rr, stress_sita;
					double sum = stress_xx[i][j][k][l] + stress_yy[i][j][k][l];
					double dif = stress_xx[i][j][k][l] - stress_yy[i][j][k][l];
					double tau2 = stress_xy[i][j][k][l] * stress_xy[i][j][k][l];
					stress_rr = sum * 0.5
								+ sqrt(dif * dif + 4 * tau2) * 0.5;
					stress_sita = sum * 0.5
								  - sqrt(dif * dif + 4 * tau2) * 0.5;					
					fprintf(fp, "% 1.10e\t% 1.10e\t% 1.10e\t% 1.10e\t% 1.10e\t% 1.10e\n", 
							calc_xi_loc[i * division_ele_xi + k], 
							calc_eta_loc[j * division_ele_eta + l],
							coord_x[i * division_ele_xi + k][j * division_ele_eta + l], 
							coord_y[i * division_ele_xi + k][j * division_ele_eta + l], 
							stress_rr, stress_sita);
				}
			}
		}
	}
	fclose(fp);

	fp = fopen("over_stress_vm_graph.txt", "a");
 	for (i = 0; i < element_n_xi; i++) {
		for (j = 0; j < element_n_eta; j++) {
			for (k = 0; k < division_ele_xi + 1; k++) {
				for (l = 0; l < division_ele_eta + 1; l++) {
					double stress_vm;
					double sum = stress_xx[i][j][k][l] + stress_yy[i][j][k][l];
					double dif = stress_xx[i][j][k][l] - stress_yy[i][j][k][l];
					double tau2 = stress_xy[i][j][k][l] * stress_xy[i][j][k][l];
					double temp1, temp2;
					temp1 = 0.5 * sum;
					temp2 = 0.5 * sqrt(dif * dif + 4 * tau2);
					stress_vm = sqrt(temp1 * temp1 + 3 * temp2 * temp2);
					fprintf(fp, "% 1.10e\t% 1.10e\t% 1.10e\t% 1.10e\t% 1.10e\n", 
							calc_xi_loc[i * division_ele_xi + k], 
							calc_eta_loc[j * division_ele_eta + l],
							coord_x[i * division_ele_xi + k][j * division_ele_eta + l], 
							coord_y[i * division_ele_xi + k][j * division_ele_eta + l], 
							stress_vm);
				}
			}
		}
	}
	fclose(fp);
}

int main(int argc, char *argv[]) {
	int i, j;

	if (argc != 8) {    //後で書き直し
		printf("Error!!\n");
		printf("Usage:\n"
			   "./NURBS_input input.txt Displacement.dat "
			   "view.dat "
               "input_local.txt "
			   "overlay_view.dat "
               "div_xi div_eta\n"
			   "\n");
		exit(1);
	}
	fp = fopen(argv[1], "r");
	if (fp == NULL) {
		printf("Error!!\n");
		printf("%s file not open!\n", argv[1]);
		exit(1);
	}
	if (argc == 6) {
		fp = fopen(argv[2], "r");
		if (fp == NULL) {
			printf("Error!!\n");
			printf("%s file not open!\n", argv[2]);
			exit(1);
		}
	}
	fclose(fp);


	division_ele_xi  = atoi(argv[argc - 2]);
	division_ele_eta = atoi(argv[argc - 1]);

	if (division_ele_xi > MAX_DIVISION) {
		printf("Error!!\n");
		printf("Too many Divsion at xi!\n"
			   "Maximum of division is %d (Now %d)\n"
			   "\n", MAX_DIVISION, division_ele_xi);
		exit(1);
	}
	if (division_ele_eta > MAX_DIVISION) {
		printf("Error!!\n");
		printf("Too many Divsion at eta!\n"
			   "Maximum of division is %d (Now %d)\n"
			   "\n", MAX_DIVISION, division_ele_eta);
		exit(1);
	}

	//ローカルメッシュの情報取得
	GetLocData (argv[4]);
   
    int patch_n_loc, patch_n_glo;	//パッチ番号

	ReadFile (argv[1], argv[2]);
	fp = fopen(argv[argc - 5], "w");
	fprintf(fp, "%d\t%d\t%d\n",
			fields_flag, division_ele_xi, division_ele_eta);
	fclose(fp);

	n_patch_glo = patch_n - n_patch_loc;

    //for s-IGA
	//重ね合わせ結果出力のためのoverlay_view.dat作成
	fp = fopen(argv[argc - 3], "w");
	fprintf(fp, "%d\t%d\t%d\n",
				fields_flag, division_ele_xi, division_ele_eta);
	fclose(fp);

	//グラフ作成のための出力
	fp = fopen("disp_graph.txt", "w");
	fprintf(fp, "patch_n\tx\ty\tdisp_x\tdisp_y\n");
	fclose(fp);

	fp = fopen("stress_y_graph.txt", "w");
	fprintf(fp, "patch_n\tx\ty\tstress_yy\n");
	fclose(fp);

	fp = fopen("stress_vm_graph.txt", "w");
	fprintf(fp, "xi\teta\tx\ty\tstress_vm\n");
	fclose(fp);

	fp = fopen("over_disp_graph.txt", "w");
	fprintf(fp, "patch_n\tx\ty\tdisp_x\tdisp_y\n");
	fclose(fp);

	fp = fopen("over_stress_x_graph.txt", "w");
	fprintf(fp, "x\ty\tstress_xx\n");
	fclose(fp);

	fp = fopen("over_stress_y_graph.txt", "w");
	fprintf(fp, "x\ty\tstress_yy\n");
	fclose(fp);

	fp = fopen("over_stress_r_graph.txt", "w");
	fprintf(fp, "xi\teta\tx\ty\tstress_rr\tstress_sita\n");
	fclose(fp);

	fp = fopen("over_stress_vm_graph.txt", "w");
	fprintf(fp, "xi\teta\tx\ty\tstress_vm\n");
	fclose(fp);

	for (i = 0; i < patch_n; i++) {

		fp = fopen("stress_vm_graph.txt", "a");
		fprintf(fp, "\npatch_n;%d\n\n",i);
		fclose(fp);	

		graph_patch_n = i;

		printf("----------Start calculation at patch %d----------\n\n", i);
		Calculation(argv[argc - 5],
					order_xi[i], order_eta[i],
					knot_n_xi[i], knot_n_eta[i],
					cntl_p_n_xi[i], cntl_p_n_eta[i],
					knot_vec_xi[i], knot_vec_eta[i],
					cntl_px[i], cntl_py[i],
					disp_cntl_px[i], disp_cntl_py[i],
					weight[i]);
		printf("-----------End calculation at patch %d-----------\n\n", i);
	
		if (i >= n_patch_glo)	//ローカル上のパッチに対しては重合計算行う
		{
			patch_n_loc = i;
	        printf("----------Start overlay calculation at patch %d in LOCAL patch----------\n\n", i);
	        for (j = 0; j < n_patch_glo; j++)
    	    {
				patch_n_glo = j;
				//printf("patch_n_loc: %d \tpatch_n_glo: %d\n",
				//		  patch_n_loc, patch_n_glo);
				Calculation_overlay(argv[argc - 3],
	                                order_xi[patch_n_loc],order_eta[patch_n_loc],
					            	knot_n_xi[patch_n_loc], knot_n_eta[patch_n_loc],
					            	cntl_p_n_xi[patch_n_loc], cntl_p_n_eta[patch_n_loc],
					            	knot_vec_xi[patch_n_loc], knot_vec_eta[patch_n_loc],
					            	cntl_px[patch_n_loc], cntl_py[patch_n_loc],
					            	disp_cntl_px[patch_n_loc], disp_cntl_py[patch_n_loc],
					            	weight[patch_n_loc],
   	                            	order_xi[patch_n_glo],order_eta[patch_n_glo],
					            	knot_n_xi[patch_n_glo], knot_n_eta[patch_n_glo],
					            	cntl_p_n_xi[patch_n_glo], cntl_p_n_eta[patch_n_glo],
					            	knot_vec_xi[patch_n_glo], knot_vec_eta[patch_n_glo],
					            	cntl_px[patch_n_glo], cntl_py[patch_n_glo],
					            	disp_cntl_px[patch_n_glo], disp_cntl_py[patch_n_glo],
				    	        	weight[patch_n_glo]);

        	}
		}	
	}
    
	return 0;
}

/*
	fp = fopen("Virtual_Crack_Extension.dat", "w");
	fprintf(fp, "label=Virtual_Crack_Extension\n");
	fprintf(fp, "num_items=%d\n", Total_Control_Point);
	fprintf( fp, "\n");
	for (j = 0; j < Total_Control_Point; j++ ) {
		fprintf(fp, "%d:	%le %le ", j,
				Virtual_Crack_Extension_Ct_Pt[j][0],
				Virtual_Crack_Extension_Ct_Pt[j][1]);
		fprintf(fp, "\n");
	}
	fclose(fp);
	*/