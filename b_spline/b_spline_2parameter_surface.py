import numpy as np
import matplotlib.pyplot as plt
import math
import function_of_NURBS as fn

# Define color vector
color = np.array(["r", "g", "b", "c", "m", "y", "k"])

# Define control points surface (for i (for j))
CP_matrix = np.array([[ 0.,   0. ],
                      [-1.,   0. ],
                      [-2.,   0. ],
                      [ 0.,   1. ],
                      [-1.,   2. ],
                      [-2.,   2. ],
                      [ 1.,   1.5],
                      [ 1.,   4. ],
                      [ 1.,   5. ],
                      [ 3.,   1.5],
                      [ 3.,   4. ],
                      [ 3.,   5. ]])

# Define 刻み幅
delta = np.array([20, 10])

# Define polynomial order:n
n = np.array([2, 2])   # n次のB-スプライン曲線

# (ξ, η)方向のコントロールポイントの数
l_i = 4
l_j = 3
l = np.array([l_i, l_j])

# reshape CP_surface
CP_surface = np.reshape(CP_matrix, (l_i, l_j, 2))

# Define number of knots 各方向ノットの個数
m = np.array([l_i+n[0]+1, l_j+n[1]+1])

# Difine knot vector
knot_i = fn.def_knot(m[0], n[0])
knot_j = fn.def_knot(m[1], n[1])

# 変数宣言
N = np.zeros((n[0]+1, delta[0], l_i))
M = np.zeros((n[1]+1, delta[1], l_j))

# 基底関数の計算
N = fn.basisfunction_return_N(N, delta[0], knot_i, l_i, m[0], n[0])
M = fn.basisfunction_return_N(M, delta[1], knot_j, l_j, m[1], n[1])

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
        for p in range(l_i):
            for q in range(l_j):
                Sx += N[n[0]][i][p] * M[n[1]][j][q] * CP_surface[p][q][0]
                Sy += N[n[0]][i][p] * M[n[1]][j][q] * CP_surface[p][q][1]
        Sx_vec[i][j] = Sx
        Sy_vec[i][j] = Sy
ax1.plot(Sx_vec, Sy_vec, c=color[0], marker="", linewidth=0.7)
ax1.plot(Sx_vec.T, Sy_vec.T, c=color[0], marker="", linewidth=0.7)

# コントロールポイントの描写
ax1.scatter(CP_matrix[:, 0], CP_matrix[:, 1], c=color[2], s=10)
x = np.zeros((l_i, l_j))
y = np.zeros((l_i, l_j))
for i in range(l_j):
    x[:, i] = CP_surface[:, i, 0]
    y[:, i] = CP_surface[:, i, 1]
ax1.plot(x, y, c=color[2], marker="", linewidth=1)
ax1.plot(x.T, y.T, c=color[2], marker="", linewidth=1)

# 描写
ax1.set_aspect('equal', adjustable='box')
ax1.set_axisbelow(True)
fig.set_figheight(9)
fig.set_figwidth(12)
ax1.grid()
ax1.set_xlim(-3, 4)
ax1.set_ylim(-1, 6)
plt.show()
