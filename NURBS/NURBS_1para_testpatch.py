import numpy as np
import matplotlib.pyplot as plt
import math
import function_of_NURBS as fn

# Define color vector
color = np.array(["r", "g", "b", "c", "m", "y", "k"])

# Define control points:CP
CP = np.array([[1.0, 0.0],
               [1.0, 1.0],
               [0.0, 1.0]])

# Define polynomial order
n = 2
l = CP.shape[0]  # 制御点の個数
m = l + n + 1   # ノットの個数

# Define knot vector
make_C0_CP = np.array([])  # C0連続にするコントロールポイント番号1個のみ(2個以上はバグる)，かつn=2のみ使える
knot = fn.def_knot_C0(m, n, make_C0_CP)

# knot insertion B
insert_knot = np.array([])  # 挿入するノットの値(0.0 < insert_knot < 1.0)
CP, l, m, knot = fn.knot_insertion_B(CP, n, l, knot, insert_knot)

# Define 刻み幅
delta = 1000

# 変数宣言
N = np.zeros((n+1, delta, l))
Cx_vec = np.zeros((delta))
Cy_vec = np.zeros((delta))
x_vec = np.zeros((delta, l))
y_vec = np.zeros((delta, l))
x_bool_vec = np.zeros((delta, l), dtype=bool)
y_bool_vec = np.zeros((delta, l), dtype=bool)

# 描写
fig = plt.figure()
# ax1 = fig.add_subplot(2, 1, 1)
ax2 = fig.add_subplot(1, 1, 1)

# 基底関数の計算
N = fn.basisfunction_return_N(N, delta, knot, l, m, n)

# Bスプラインの描写
for i in range(delta):
    Cx = 0
    Cy = 0
    for j in range(l):
        Cx += N[n][i][j] * CP[j][0]
        Cy += N[n][i][j] * CP[j][1]
    Cx_vec[i] = Cx
    Cy_vec[i] = Cy
ax2.plot(Cx_vec, Cy_vec, c=color[0], marker="", linewidth=1.0)

# コントロールポイントの描写
for i in range(l):
    x = CP[i][0]
    y = CP[i][1]
    ax2.scatter(x, y, c=color[2], s=8)
ax2.plot(CP[:, 0], CP[:, 1], c=color[2], marker="", linewidth=1.0)

# 描写
ax2.grid()
ax2.set_xlim(0, 1.2)
ax2.set_ylim(0, 1.2)
ax2.set_aspect('equal', adjustable='box')

# Define control points
# weight value
wv = 1./math.sqrt(2)

CP_kinji_0 = np.array([[1., 0., 0., 1.],
                       [1., 1., 0., 1.*wv],
                       [0., 1., 0., 1.]])

CP_kinji = CP_kinji_0[:, :-1]
w = CP_kinji_0[:, 3:]

knot = np.array([0., 0., 0., 1., 1., 1.])

# 変数宣言
delta = 500
n = 2
l = 3
m = n + l + 1
N = np.zeros((n+1, delta, l))
R = np.zeros((delta, l))
Cx_vec_1 = np.zeros((delta))
Cy_vec_1 = np.zeros((delta))
Cy_vec_2 = np.zeros((delta))
N = fn.basisfunction_return_N(N, delta, knot, l, m, n)
R = fn.weight_basisfunction_1parameter_return_R(R, N, w, delta, n, l)
# 重み付き近似スプライン
for i in range(delta):
    Cx = 0
    Cy = 0
    for j in range(l):
        Cx += R[i][j] * CP_kinji[j][0]
        Cy += R[i][j] * CP_kinji[j][1]
    Cx_vec_1[i] = Cx
    Cy_vec_1[i] = Cy
    Cy_vec_2[i] = 5*math.sin(math.acos(Cx/5.))
ax2.plot(Cx_vec_1, Cy_vec_1, c=color[3], marker="", linewidth=1.0)

ax2.set_axisbelow(True)
ax2.set_xlabel("x")
ax2.set_ylabel("y")
fig.set_figheight(4)
fig.set_figwidth(4)
plt.show()
