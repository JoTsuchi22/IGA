import numpy as np
import matplotlib.pyplot as plt
import math
import function_of_NURBS_alt as fn

# Define color vector
color = np.array(["r", "g", "b", "c", "m", "y", "k"])

# Define control points
# weight value
th = math.pi/4.
v0 = math.cos(th)
v1 = math.sin(th)
wv = math.cos(th/2.)

CP_matrix_weight = np.array([[1.,  0., 1.],
                             [1., math.sin(th/2.), wv],
                             [v0,  v1, 1.],
                             [3.,  0., 1.],
                             [3., 1.5, 1.],
                             [3.,  3., 1.]])

# Define polynomial order:n
n = np.array([1, 2])   # n次のB-スプライン曲線

# (ξ, η)方向のコントロールポイントの数
l_i = 2
l_j = 3
l = np.array([l_i, l_j])

# Define number of knots 各方向ノットの個数
m = np.array([l[0]+n[0]+1, l[1]+n[1]+1])

# Difine knot vector
knot_i = fn.def_knot(m[0], n[0])
# knot_j = fn.def_knot(m[1], n[1])
# ノットの置き方 特殊
# knot_i = np.array([0, 0, 0, 1, 1, 1])
knot_j = np.array([0., 0., 0., 1., 1., 1.])

# CPreshape
CP_2d1w = np.reshape(CP_matrix_weight, (l[0], l[1], 3))

# オーダーエレベーション
elevation_parameter_axis = 0
elevation_degree = 1
CP_2d1w, l, m, n, knot_i, knot_j = fn.order_elevation_2p2d1w(
    CP_2d1w, l, m, n, knot_i, knot_j, elevation_degree, elevation_parameter_axis)

# オーダーエレベーション
elevation_parameter_axis = 1
elevation_degree = 0
CP_2d1w, l, m, n, knot_i, knot_j = fn.order_elevation_2p2d1w(
    CP_2d1w, l, m, n, knot_i, knot_j, elevation_degree, elevation_parameter_axis)

# autoノットインサーション
insert_parameter_axis = 0
number_of_auto_insertion = 2
CP_2d1w, l, m, knot_i, knot_j = fn.knot_insertion_C_2p2d1w(CP_2d1w, n, l, m, knot_i, knot_j,
                                                           insert_parameter_axis, number_of_auto_insertion)

# autoノットインサーション
insert_parameter_axis = 1
number_of_auto_insertion = 2
CP_2d1w, l, m, knot_i, knot_j = fn.knot_insertion_C_2p2d1w(CP_2d1w, n, l, m, knot_i, knot_j,
                                                           insert_parameter_axis, number_of_auto_insertion)

# ノットインサーションB
insert_parameter_axis = 1
insert_knot = np.array([])
CP_2d1w, l, m, knot_i, knot_j = fn.knot_insertion_B_2p2d1w(CP_2d1w, n, l, m, knot_i, knot_j,
                                                           insert_parameter_axis, insert_knot)

# コントロールポイントreshape
CP_matrix = np.reshape(
    CP_2d1w, [int((CP_2d1w.shape[0]*CP_2d1w.shape[1]*CP_2d1w.shape[2])/3.), 3])[:, :-1]
weight = np.reshape(
    CP_2d1w, [int((CP_2d1w.shape[0]*CP_2d1w.shape[1]*CP_2d1w.shape[2])/3.), 3])[:, 2:]

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

CP_matrix = fn.affine_transformation_2D(CP_matrix, CP_matrix.shape[0], stretch_x, stretch_y, stretch_z,
                                        trans_x, trans_y, trans_z, theta_x, theta_y, theta_z, shear_x, shear_y)

# reshape CP_2D
CP_2D = np.reshape(CP_matrix, (l[0], l[1], 2))

# reshape weight
w = np.reshape(weight, (l[0], l[1]))

# Define 刻み幅
delta = np.array([l[0]-1, l[1]-1])
# delta = np.array([m[0], m[1]])
# delta = np.array([6, 6])

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

# 描写
ax1.set_aspect('equal', adjustable='box')
ax1.set_axisbelow(True)
fig.set_figheight(9)
fig.set_figwidth(12)
ax1.grid()
ax1.set_xlim(0, 5)
ax1.set_ylim(-1, 4)
plt.show()

# solid_name = "example"
# fn.make_stl_3D(solid_name, delta, Sx_vec , Sy_vec, Sz_vec)
