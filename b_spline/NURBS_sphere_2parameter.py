import numpy as np
import matplotlib.pyplot as plt
import math
import function_of_NURBS as fn

# Define color vector
color = np.array(["r", "g", "b", "c", "m", "y", "k"])

# Define control points
# weight value
wv = 1 / math.sqrt(2)

CP_matrix_weight = np.array([[0., 0., 0., 1.],
                             [0., 0., 0., wv],
                             [0., 0., 0., 1.],
                             [0., 0., 0., wv],
                             [0., 0., 0., 1.],
                             [0., 0., 0., wv],
                             [0., 0., 0., 1.],
                             [0., 0., 0., wv],
                             [0., 0., 0., 1.],
                             [0.,  1.,  0.,  1.*wv],
                             [0.,  1.,  1.,  wv*wv],
                             [0.,  0.,  1.,  1.*wv],
                             [0., -1.,  1.,  wv*wv],
                             [0., -1.,  0.,  1.*wv],
                             [0., -1., -1.,  wv*wv],
                             [0.,  0., -1.,  1.*wv],
                             [0.,  1., -1.,  wv*wv],
                             [0.,  1.,  0.,  1.*wv],
                             [1.,  1.,  0.,  1.],
                             [1.,  1.,  1.,  wv],
                             [1.,  0.,  1.,  1.],
                             [1., -1.,  1.,  wv],
                             [1., -1.,  0.,  1.],
                             [1., -1., -1.,  wv],
                             [1.,  0., -1.,  1.],
                             [1.,  1., -1.,  wv],
                             [1.,  1.,  0.,  1.],
                             [2.,  1.,  0.,  1.*wv],
                             [2.,  1.,  1.,  wv*wv],
                             [2.,  0.,  1.,  1.*wv],
                             [2., -1.,  1.,  wv*wv],
                             [2., -1.,  0.,  1.*wv],
                             [2., -1., -1.,  wv*wv],
                             [2.,  0., -1.,  1.*wv],
                             [2.,  1., -1.,  wv*wv],
                             [2.,  1.,  0.,  1.*wv],
                             [2., 0., 0., 1.],
                             [2., 0., 0., wv],
                             [2., 0., 0., 1.],
                             [2., 0., 0., wv],
                             [2., 0., 0., 1.],
                             [2., 0., 0., wv],
                             [2., 0., 0., 1.],
                             [2., 0., 0., wv],
                             [2., 0., 0., 1.]])

CP_matrix = CP_matrix_weight[:, :-1]
weight = CP_matrix_weight[:, 3:]

# affine transformation (for CP)
stretch_x = 2.0
stretch_y = 2.0
stretch_z = 2.0

trans_x = 0.0
trans_y = 0.0
trans_z = 0.0

theta_x = 0.0
theta_y = 0.0
theta_z = 0.0

shear_x = 0.0
shear_y = 0.0

CP_matrix = fn.affine_transformation_3D(CP_matrix, CP_matrix.shape[0], stretch_x, stretch_y, stretch_z,
                                         trans_x, trans_y, trans_z, theta_x, theta_y, theta_z, shear_x, shear_y)

# Define 刻み幅
delta = np.array([60, 60])

# Define polynomial order:n
n = np.array([2, 2])   # n次のB-スプライン曲線

# (ξ, η)方向のコントロールポイントの数
l_i = 5
l_j = 9
l = np.array([l_i, l_j])

# reshape CP_3D
CP_3D = np.reshape(CP_matrix, (l_i, l_j, 3))

# reshape weight
w = np.reshape(weight, (l_i, l_j))

# Define number of knots 各方向ノットの個数
m = np.array([l_i+n[0]+1, l_j+n[1]+1])

# Difine knot vector
# knot_i = fn.def_knot(m[0], n[0])
# knot_j = fn.def_knot(m[1], n[1])
# ノットの置き方 特殊
knot_i = np.array([0, 0, 0, 1, 1, 2, 2, 2])
knot_j = np.array([0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4])

# 変数宣言
N = np.zeros((n[0]+1, delta[0], l_i))
M = np.zeros((n[1]+1, delta[1], l_j))
R = np.zeros((delta[0], delta[1], l_i, l_j))

# 基底関数の計算
N = fn.basisfunction_return_N(N, delta[0], knot_i, l_i, m[0], n[0])
M = fn.basisfunction_return_N(M, delta[1], knot_j, l_j, m[1], n[1])

# 重み付き基底関数の計算
R = fn.weight_basisfunction_2parameter_return_R(R, N, M, w, delta, n, l)

# 描写
fig = plt.figure(figsize=(10, 10))
ax1 = fig.add_subplot(111, projection='3d')

# Bスプラインの描写
Sx_vec = np.zeros((delta[0], delta[1]))
Sy_vec = np.zeros((delta[0], delta[1]))
Sz_vec = np.zeros((delta[0], delta[1]))
for i in range(delta[0]):
    for j in range(delta[1]):
        Sx = 0
        Sy = 0
        Sz = 0
        for p in range(l_i):
            for q in range(l_j):
                Sx += R[i][j][p][q] * CP_3D[p][q][0]
                Sy += R[i][j][p][q] * CP_3D[p][q][1]
                Sz += R[i][j][p][q] * CP_3D[p][q][2]
        Sx_vec[i][j] = Sx
        Sy_vec[i][j] = Sy
        Sz_vec[i][j] = Sz

# ワイヤーフレーム表示
# for i in range(delta[0]):
#     ax1.plot(Sx_vec[i,:], Sy_vec[i,:], Sz_vec[i,:], c=color[0], linewidth=0.3)
# for j in range(delta[1]):
#     ax1.plot(Sx_vec[:,j], Sy_vec[:,j], Sz_vec[:,j], c=color[0], linewidth=0.3)

# メッシュ表示
ax1.plot_surface(Sx_vec[:, :], Sy_vec[:, :],Sz_vec[:, :], cmap="viridis", alpha=0.5)

# # # コントロールポイントの描写
x = np.zeros((l_i, l_j))
y = np.zeros((l_i, l_j))
z = np.zeros((l_i, l_j))
for i in range(l_i):
    for j in range(l_j):
        x[i, j] = CP_3D[i, j, 0]
        y[i, j] = CP_3D[i, j, 1]
        z[i, j] = CP_3D[i, j, 2]
for i in range(l_i):
    ax1.plot(x[i,:], y[i,:], z[i,:], c=color[2], linewidth=1)
for j in range(l_j):
    ax1.plot(x[:,j], y[:,j], z[:,j], c=color[2], linewidth=1)

# 描写
ax1.set_box_aspect((1, 1, 1))
ax1.set_axisbelow(True)
ax1.grid()
ax1.set_xlim(-5, 5)
ax1.set_ylim(-5, 5)
ax1.set_zlim(-5, 5)
fig.set_figheight(9)
fig.set_figwidth(12)
plt.show()

solid_name = "sphere_2parameter"
fn.make_stl_2D(solid_name, delta, Sx_vec , Sy_vec, Sz_vec)