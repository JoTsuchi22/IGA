import numpy as np
import matplotlib.pyplot as plt
import math
import b_spline_function as bpf

# Define color vector
color = np.array(["r", "g", "b", "c", "m", "y", "k"])

# Define control points
# weight value
wv = 1 / math.sqrt(2)

CP_matrix_weight = np.array([[ 1.,   0.,   4. , 1.],
                             [ 1.,   0.,   2.5, 1.],
                             [ 1.,   0.,   1. , 1.],
                             [ 1.,   0.,  -1. , wv],
                             [ 3.,   0.,  -1. , 1.],
                             [ 1.,   1.,   4. , 1.*wv],
                             [ 1.,   1.,   2.5, 1.*wv],
                             [ 1.,   1.,   1. , 1.*wv],
                             [ 1.,   1.,  -1. , wv*wv],
                             [ 3.,   1.,  -1. , 1.*wv],
                             [ 0.,   1.,   4. , 1.],
                             [ 0.,   1.,   2.5, 1.],
                             [ 0.,   1.,   1. , 1.],
                             [ 0.,   1.,  -2. , wv],
                             [ 3.,   1.,  -2. , 1.],
                             [-1.,   1.,   4. , 1.*wv],
                             [-1.,   1.,   2.5, 1.*wv],
                             [-1.,   1.,   1. , 1.*wv],
                             [-1.,   1.,  -3. , wv*wv],
                             [ 3.,   1.,  -3. , 1.*wv],
                             [-1.,   0.,   4. , 1.],
                             [-1.,   0.,   2.5, 1.],
                             [-1.,   0.,   1. , 1.],
                             [-1.,   0.,  -3. , wv],
                             [ 3.,   0.,  -3. , 1.],
                             [-1.,  -1.,   4. , 1.*wv],
                             [-1.,  -1.,   2.5, 1.*wv],
                             [-1.,  -1.,   1. , 1.*wv],
                             [-1.,  -1.,  -3. , wv*wv],
                             [ 3.,  -1.,  -3. , 1.*wv],
                             [ 0.,  -1.,   4. , 1.],
                             [ 0.,  -1.,   2.5, 1.],
                             [ 0.,  -1.,   1. , 1.],
                             [ 0.,  -1.,  -2. , wv],
                             [ 3.,  -1.,  -2. , 1.],
                             [ 1.,  -1.,   4. , 1.*wv],
                             [ 1.,  -1.,   2.5, 1.*wv],
                             [ 1.,  -1.,   1. , 1.*wv],
                             [ 1.,  -1.,  -1. , wv*wv],
                             [ 3.,  -1.,  -1. , 1.*wv],
                             [ 1.,   0.,   4. , 1.],
                             [ 1.,   0.,   2.5, 1.],
                             [ 1.,   0.,   1. , 1.],
                             [ 1.,   0.,  -1. , wv],
                             [ 3.,   0.,  -1. , 1.],
                             [ 2.,   0.,   4. , 1.],
                             [ 2.,   0.,   2.5, 1.],
                             [ 2.,   0.,   1. , 1.],
                             [ 2.,   0.,   0. , wv],
                             [ 3.,   0.,   0. , 1.],
                             [ 2.,   2.,   4. , 1.*wv],
                             [ 2.,   2.,   2.5, 1.*wv],
                             [ 2.,   2.,   1. , 1.*wv],
                             [ 2.,   2.,   0. , wv*wv],
                             [ 3.,   2.,   0. , 1.*wv],
                             [ 0.,   2.,   4. , 1.],
                             [ 0.,   2.,   2.5, 1.],
                             [ 0.,   2.,   1. , 1.],
                             [ 0.,   2.,  -2. , wv],
                             [ 3.,   2.,  -2. , 1.],
                             [-2.,   2.,   4. , 1.*wv],
                             [-2.,   2.,   2.5, 1.*wv],
                             [-2.,   2.,   1. , 1.*wv],
                             [-2.,   2.,  -4. , wv*wv],
                             [ 3.,   2.,  -4. , 1.*wv],
                             [-2.,   0.,   4. , 1.],
                             [-2.,   0.,   2.5, 1.],
                             [-2.,   0.,   1. , 1.],
                             [-2.,   0.,  -4. , wv],
                             [ 3.,   0.,  -4. , 1.],
                             [-2.,  -2.,   4. , 1.*wv],
                             [-2.,  -2.,   2.5, 1.*wv],
                             [-2.,  -2.,   1. , 1.*wv],
                             [-2.,  -2.,  -4. , wv*wv],
                             [ 3.,  -2.,  -4. , 1.*wv],
                             [ 0.,  -2.,   4. , 1.],
                             [ 0.,  -2.,   2.5, 1.],
                             [ 0.,  -2.,   1. , 1.],
                             [ 0.,  -2.,  -2. , wv],
                             [ 3.,  -2.,  -2. , 1.],
                             [ 2.,  -2.,   4. , 1.*wv],
                             [ 2.,  -2.,   2.5, 1.*wv],
                             [ 2.,  -2.,   1. , 1.*wv],
                             [ 2.,  -2.,   0. , wv*wv],
                             [ 3.,  -2.,   0. , 1.*wv],
                             [ 2.,   0.,   4. , 1.],
                             [ 2.,   0.,   2.5, 1.],
                             [ 2.,   0.,   1. , 1.],
                             [ 2.,   0.,   0. , wv],
                             [ 3.,   0.,   0. , 1.]])

CP_matrix = CP_matrix_weight[:, :-1]
weight = CP_matrix_weight[:, 3:]

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

CP_matrix = bpf.affine_transformation_3D(CP_matrix, CP_matrix.shape[0], stretch_x, stretch_y, stretch_z,
                                         trans_x, trans_y, trans_z, theta_x, theta_y, theta_z, shear_x, shear_y)

# Define 刻み幅
delta = np.array([6, 25, 25])

# Define polynomial order:n
n = np.array([1, 2, 2])   # n次のB-スプライン曲線

# (ξ, η)方向のコントロールポイントの数
l_i = 2
l_j = 9
l_k = 5
l = np.array([l_i, l_j, l_k])

# reshape CP_3D
CP_3D = np.reshape(CP_matrix, (l_i, l_j, l_k, 3))

# reshape weight
w = np.reshape(weight, (l_i, l_j, l_k))

# Define number of knots 各方向ノットの個数
m = np.array([l_i+n[0]+1, l_j+n[1]+1, l_k+n[2]+1])

# Difine knot vector
# knot_i = bpf.def_knot(m[0], n[0])
# knot_j = bpf.def_knot(m[1], n[1])
# knot_k = bpf.def_knot(m[2], n[2])
# ノットの置き方 特殊
knot_i = np.array([0, 0, 1, 1])
knot_j = np.array([0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4])
knot_k = np.array([0, 0, 0, 1, 1, 2, 2, 2])

# 変数宣言
N = np.zeros((n[0]+1, delta[0], l_i))
M = np.zeros((n[1]+1, delta[1], l_j))
L = np.zeros((n[2]+1, delta[2], l_k))
R = np.zeros((delta[0], delta[1], delta[2], l_i, l_j, l_k))

# 基底関数の計算
N = bpf.basisfunction_return_N(N, delta[0], knot_i, l_i, m[0], n[0])
M = bpf.basisfunction_return_N(M, delta[1], knot_j, l_j, m[1], n[1])
L = bpf.basisfunction_return_N(L, delta[2], knot_k, l_k, m[2], n[2])

# 重み付き基底関数の計算
R = bpf.weight_basisfunction_3parameter_return_R(R, N, M, L, w, delta, n, l)

# 描写
fig = plt.figure(figsize=(10, 10))
ax1 = fig.add_subplot(111, projection='3d')

# Bスプラインの描写
Sx_vec = np.zeros((delta[0], delta[1], delta[2]))
Sy_vec = np.zeros((delta[0], delta[1], delta[2]))
Sz_vec = np.zeros((delta[0], delta[1], delta[2]))
for i in range(delta[0]):
    for j in range(delta[1]):
        for k in range(delta[2]):
            Sx = 0
            Sy = 0
            Sz = 0
            for p in range(l_i):
                for q in range(l_j):
                    for r in range(l_k):
                        Sx += R[i][j][k][p][q][r] * CP_3D[p][q][r][0]
                        Sy += R[i][j][k][p][q][r] * CP_3D[p][q][r][1]
                        Sz += R[i][j][k][p][q][r] * CP_3D[p][q][r][2]
            Sx_vec[i][j][k] = Sx
            Sy_vec[i][j][k] = Sy
            Sz_vec[i][j][k] = Sz

# ワイヤーフレーム表示
# for i in range(delta[0]):
#     for j in range(delta[1]):
#         ax1.plot(Sx_vec[i, j, :], Sy_vec[i, j, :],
#                  Sz_vec[i, j, :], c=color[0], linewidth=0.3)
# for i in range(delta[0]):
#     for k in range(delta[2]):
#         ax1.plot(Sx_vec[i, :, k], Sy_vec[i, :, k],
#                  Sz_vec[i, :, k], c=color[0], linewidth=0.3)
# for j in range(delta[1]):
#     for k in range(delta[2]):
#         ax1.plot(Sx_vec[:, j, k], Sy_vec[:, j, k],
#                  Sz_vec[:, j, k], c=color[0], linewidth=0.3)

# メッシュ表示
ax1.plot_surface(Sx_vec[:, :, 0], Sy_vec[:, :, 0],
                 Sz_vec[:, :, 0], cmap="viridis", alpha=0.5)
ax1.plot_surface(Sx_vec[:, :, -1], Sy_vec[:, :, -1],
                 Sz_vec[:, :, -1], cmap="viridis", alpha=0.5)
ax1.plot_surface(Sx_vec[0, :, :], Sy_vec[0, :, :],
                 Sz_vec[0, :, :], cmap="viridis", alpha=0.5)
ax1.plot_surface(Sx_vec[-1, :, :], Sy_vec[-1, :, :],
                 Sz_vec[-1, :, :], cmap="viridis", alpha=0.5)
# ax1.plot_surface(Sx_vec[:, 0, :], Sy_vec[:, 0, :],
#                  Sz_vec[:, 0, :], cmap="viridis", alpha=0.5)
# ax1.plot_surface(Sx_vec[:, -1, :], Sy_vec[:, -1, :],
#                  Sz_vec[:, -1, :], cmap="viridis", alpha=0.5)


# # # コントロールポイントの描写
x = np.zeros((l_i, l_j, l_k))
y = np.zeros((l_i, l_j, l_k))
z = np.zeros((l_i, l_j, l_k))
for i in range(l_i):
    for j in range(l_j):
        for k in range(l_k):
            x[i, j, k] = CP_3D[i, j, k, 0]
            y[i, j, k] = CP_3D[i, j, k, 1]
            z[i, j, k] = CP_3D[i, j, k, 2]
for i in range(l_i):
    for j in range(l_j):
        ax1.plot(x[i, j, :], y[i, j, :], z[i, j, :], c=color[2], linewidth=0.6)
for i in range(l_i):
    for k in range(l_k):
        ax1.plot(x[i, :, k], y[i, :, k], z[i, :, k], c=color[2], linewidth=0.6)
for j in range(l_j):
    for k in range(l_k):
        ax1.plot(x[:, j, k], y[:, j, k], z[:, j, k], c=color[2], linewidth=0.6)

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

solid_name = "example"
bpf.make_stl_3D(solid_name, delta, Sx_vec , Sy_vec, Sz_vec)