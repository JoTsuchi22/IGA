import numpy as np
import matplotlib.pyplot as plt
import math
import function_of_NURBS as fn

# Define color vector
color = np.array(["r", "g", "b", "c", "m", "y", "k"])

# Define control points:CP
CP = np.array([[1.0, 2.0],
               [2.0, 1.0],
               [3.0, 1.0],
               [3.0, 3.0],
               [5.0, 3.0],
               [2.0, 5.0]])

# CP = np.array([[1., 2.],
#                [6., 2.]])

# CP = np.array([[0., 0.],
#                [0., 1.],
#                [1., 1.],
#                [2., 1.],
#                [2., 0.],
#                [2.,-1.],
#                [1.,-1.],
#                [0.,-1.],
#                [0., 0.]])

# Define polynomial order
n = 2
l = CP.shape[0]  # 制御点の個数
m = l + n + 1   # ノットの個数

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

CP = fn.affine_transformation_2D(CP, l, stretch_x, stretch_y, stretch_z,
                                  trans_x, trans_y, trans_z, theta_x, theta_y, theta_z, shear_x, shear_y)

# Define knot vector
make_C0_CP = np.array([])  # C0連続にするコントロールポイント番号1個のみ(2個以上はバグる)，かつn=2のみ使える
knot = fn.def_knot_C0(m, n, make_C0_CP)
# knot = np.array([0.,0.,0.,0.25,0.25,0.5,0.5,0.75,0.75,1.,1.,1.])

# knot removal
removal_knot = np.array([])  # 除去するノットの値(0.0 < insert_knot < 1.0)
CP, l, m, knot = fn.knot_removal(CP, n, l, knot, removal_knot)

# order elevation
elevation_degree = 1  # order elevationを行う回数
CP, l, m, n, knot = fn.order_elevation(CP, n, l, knot, elevation_degree)

# knot insertion A
new_knot_position = np.array([])  # ノットを挿入するコントロールポイント番号，0, 1, 2...
CP, l, m, knot = fn.knot_insertion_A(CP, n, l, knot, new_knot_position)

# knot insertion A
new_knot_position = np.array([])  # ノットを挿入するコントロールポイント番号，0, 1, 2...
CP, l, m, knot = fn.knot_insertion_A(CP, n, l, knot, new_knot_position)

# knot insertion B
insert_knot = np.array([0.25, 0.5, 0.75])  # 挿入するノットの値(0.0 < insert_knot < 1.0)
CP, l, m, knot = fn.knot_insertion_B(CP, n, l, knot, insert_knot)

# knot insertion C
number_of_auto_insertion = 0    # autoノットインサーションの回数
CP, l, m, knot = fn.knot_insertion_C(CP, n, l, m, knot, number_of_auto_insertion)

# print(CP)
print("CP = ")
print(CP)
print("knot = ")
print(knot)

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
# # ax2 = fig.add_subplot(2, 2, 3)    # ax4なしの時
# # ax3 = fig.add_subplot(2, 2, 4)    # ax4なしの時
# ax2 = fig.add_subplot(2, 3, 4)      # ax4ありの時
# ax3 = fig.add_subplot(2, 3, 5)      # ax4ありの時
# ax4 = fig.add_subplot(2, 3, 6)      # ax4ありの時

# ax1 = fig.add_subplot(2, 1, 1)
ax2 = fig.add_subplot(1, 1, 1)

# 基底関数の計算
N = fn.basisfunction_return_N(N, delta, knot, l, m, n)

# 基底関数の描写
for i in range(delta):
    for j in range(l):
        nowknot = (i/(delta - 1)) * knot[m-1]
        x = nowknot
        y = N[n][i][j]
        x_vec[i][j] = x
        y_vec[i][j] = y
        x_bool_vec[i][j] = True
        y_bool_vec[i][j] = True
#         if i == delta-1:
#             ax1.plot(x_vec[:, j][x_bool_vec[:, j]], y_vec[:, j]
#                      [y_bool_vec[:, j]], c=color[j % 7], marker="", linewidth=1.0)
# ax1.grid()
# ax1.set_xlim(0, 1)
# ax1.set_ylim(0, 1)
# ax1.set_xlabel("ξ [-]")

# Bスプラインの描写
for i in range(delta):
    Cx = 0
    Cy = 0
    for j in range(l):
        Cx += N[n][i][j] * CP[j][0]
        Cy += N[n][i][j] * CP[j][1]
    Cx_vec[i] = Cx
    Cy_vec[i] = Cy
ax2.plot(Cx_vec, Cy_vec, c=color[0], marker="", linewidth=1)

# コントロールポイントの描写
for i in range(l):
    x = CP[i][0]
    y = CP[i][1]
    ax2.scatter(x, y, c=color[2], s=5)
ax2.plot(CP[:, 0], CP[:, 1], c=color[2], marker="", linewidth=1)

# 描写
ax2.grid()
ax2.set_xlim(0, 6)
ax2.set_ylim(0, 6)
ax2.set_aspect('equal', adjustable='box')
ax2.set_xlabel("x")
ax2.set_ylabel("y")

# ax1.set_axisbelow(True)
ax2.set_axisbelow(True)

# fig.set_figheight(9)
# fig.set_figwidth(15)
fig.set_figheight(3)
fig.set_figwidth(3)
plt.show()
