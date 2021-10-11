import numpy as np
import matplotlib.pyplot as plt
import math
import function_of_NURBS_alt as fn

# output file name
file_name = "enk_loc"

# Define color vector
color = np.array(["r", "g", "b", "c", "m", "y", "k"])

# Define control points
# weight value
pi = math.pi
wv = math.cos(pi/4.)

CP_matrix_weight = np.array([[1.,  0., 1.],
                             [2.,  0., 1.],
                             [1.,  1., wv],
                             [2.,  2., wv],
                             [0.,  1., 1.],
                             [0.,  2., 1.]])

# Define polynomial order:n
n_xi =  1
n_eta = 2

# (ξ, η)方向のコントロールポイントの数
l_xi = 2
l_eta = 3

# n, lとxi, etaの関係
n = np.array([n_eta, n_xi])
l = np.array([l_eta, l_xi])

# Define number of knots 各方向ノットの個数
m = np.array([l[0]+n[0]+1, l[1]+n[1]+1])

# Difine knot vector
knot_xi = fn.def_knot(m[1], n[1])
knot_eta = fn.def_knot(m[0], n[0])

# ξはm[1]，ηはm[0]
# knot_xi = fn.def_knot(m[1], n[1])
# knot_eta = fn.def_knot(m[0], n[0])

knot_i = knot_eta
knot_j = knot_xi

# affine transformation (for CP)
stretch_x = 1.0
stretch_y = 1.0
stretch_z = 1.0

trans_x = 0.0
trans_y = 0.0
trans_z = 0.0

theta_x = 0.0
theta_y = 0.0
theta_z = 0.0

shear_x = 0.0
shear_y = 0.0

CP_matrix_weight = fn.affine_transformation_2D(CP_matrix_weight, CP_matrix_weight.shape[0], stretch_x, stretch_y, stretch_z,
                                        trans_x, trans_y, trans_z, theta_x, theta_y, theta_z, shear_x, shear_y)

# CPreshape
CP_2d1w = np.reshape(CP_matrix_weight, (l[0], l[1], 3))

#difine axis
xi = 1
eta = 0

# オーダーエレベーション xi
elevation_parameter_axis = xi
elevation_degree = 1
CP_2d1w, l, m, n, knot_i, knot_j = fn.NURBS_order_elevation_2p2d1w(
    CP_2d1w, l, m, n, knot_i, knot_j, elevation_degree, elevation_parameter_axis)

# オーダーエレベーション eta
elevation_parameter_axis = eta
elevation_degree = 0
CP_2d1w, l, m, n, knot_i, knot_j = fn.NURBS_order_elevation_2p2d1w(
    CP_2d1w, l, m, n, knot_i, knot_j, elevation_degree, elevation_parameter_axis)

# autoノットインサーション xi
insert_parameter_axis = xi
number_of_auto_insertion = 0
CP_2d1w, l, m, knot_i, knot_j = fn.NURBS_knot_insertion_C_2p2d1w(CP_2d1w, n, l, m, knot_i, knot_j,
                                                           insert_parameter_axis, number_of_auto_insertion)

# autoノットインサーション eta
insert_parameter_axis = eta
number_of_auto_insertion = 1
CP_2d1w, l, m, knot_i, knot_j = fn.NURBS_knot_insertion_C_2p2d1w(CP_2d1w, n, l, m, knot_i, knot_j,
                                                           insert_parameter_axis, number_of_auto_insertion)

# ノットインサーションB eta
insert_parameter_axis = eta
insert_knot = np.array([])
CP_2d1w, l, m, knot_i, knot_j = fn.NURBS_knot_insertion_B_2p2d1w(CP_2d1w, n, l, m, knot_i, knot_j,
                                                           insert_parameter_axis, insert_knot)

# コントロールポイントreshape
CP_matrix = np.reshape(
    CP_2d1w, [int((CP_2d1w.shape[0]*CP_2d1w.shape[1]*CP_2d1w.shape[2])/3.), 3])[:, :-1]
weight = np.reshape(
    CP_2d1w, [int((CP_2d1w.shape[0]*CP_2d1w.shape[1]*CP_2d1w.shape[2])/3.), 3])[:, 2:]

# reshape CP_2D
CP_2D = np.reshape(CP_matrix, (l[0], l[1], 2))

# reshape weight
w = np.reshape(weight, (l[0], l[1]))

# Define 刻み幅
delta = np.array([l[0]-1, l[1]-1])

# 変数宣言
N = np.zeros((n[0]+1, delta[0], l[0]))
M = np.zeros((n[1]+1, delta[1], l[1]))
R = np.zeros((delta[0], delta[1], l[0], l[1]))

# 基底関数の計算
N = fn.basisfunction_return_N(N, delta[0], knot_i, l[0], m[0], n[0])
M = fn.basisfunction_return_N(M, delta[1], knot_j, l[1], m[1], n[1])

# 重み付き基底関数の計算
R = fn.weight_basisfunction_2parameter_return_R(R, N, M, w, delta, n, l)

# 描写
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)

# Bスプラインの描写
Sx_vec = np.zeros((delta[0], delta[1]))
Sy_vec = np.zeros((delta[0], delta[1]))
for i in range(delta[0]):
    for j in range(delta[1]):
        Sx = 0
        Sy = 0
        for p in range(l[0]):
            for q in range(l[1]):
                Sx += R[i][j][p][q] * CP_2D[p][q][0]
                Sy += R[i][j][p][q] * CP_2D[p][q][1]
        Sx_vec[i][j] = Sx
        Sy_vec[i][j] = Sy
ax1.plot(Sx_vec, Sy_vec, c=color[0], marker="", linewidth=0.7)
ax1.plot(Sx_vec.T, Sy_vec.T, c=color[0], marker="", linewidth=0.7)

# コントロールポイントの描写
ax1.scatter(CP_matrix[:, 0], CP_matrix[:, 1], c=color[2], s=10)
x = np.zeros((l[0], l[1]))
y = np.zeros((l[0], l[1]))
for i in range(l[1]):
    x[:, i] = CP_2D[:, i, 0]
    y[:, i] = CP_2D[:, i, 1]
ax1.plot(x, y, c=color[2], marker="", linewidth=1)
ax1.plot(x.T, y.T, c=color[2], marker="", linewidth=1)

# テキストファイル出力
fn.output_2p2d1w_txt(file_name, n, m, l, knot_i, knot_j, CP_2d1w)

# 描写
ax1.set_aspect('equal', adjustable='box')
ax1.set_axisbelow(True)
fig.set_figheight(9)
fig.set_figwidth(12)
ax1.grid()
ax1.set_xlim(-1, 4)
ax1.set_ylim(-1, 4)
plt.show()

# solid_name = "example"
# fn.make_stl_3D(solid_name, delta, Sx_vec , Sy_vec, Sz_vec)