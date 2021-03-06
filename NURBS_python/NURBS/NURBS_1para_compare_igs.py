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
ax1 = fig.add_subplot(2, 1, 1)
ax2 = fig.add_subplot(2, 2, 3)
ax3 = fig.add_subplot(2, 2, 4)

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
        if i == delta-1:
            ax1.plot(x_vec[:, j][x_bool_vec[:, j]], y_vec[:, j]
                     [y_bool_vec[:, j]], c=color[j % 7], marker="", linewidth=0.5)
ax1.grid()
ax1.set_xlim(0, 1)
ax1.set_ylim(0, 1)

# Bスプラインの描写
for i in range(delta):
    Cx = 0
    Cy = 0
    for j in range(l):
        Cx += N[n][i][j] * CP[j][0]
        Cy += N[n][i][j] * CP[j][1]
    Cx_vec[i] = Cx
    Cy_vec[i] = Cy
ax2.plot(Cx_vec, Cy_vec, c=color[0], marker="", linewidth=0.5)

# コントロールポイントの描写
for i in range(l):
    x = CP[i][0]
    y = CP[i][1]
    ax2.scatter(x, y, c=color[2], s=5)
ax2.plot(CP[:, 0], CP[:, 1], c=color[2], marker="", linewidth=0.5)

# 描写
ax2.grid()
ax2.set_xlim(0, 2)
ax2.set_ylim(0, 2)
ax2.set_aspect('equal', adjustable='box')

# 描写
ax3.plot(Cx_vec, Cy_vec, c=color[0], marker="", linewidth=0.5)
xi = np.unique(knot)
N_xi = np.zeros((n+1, xi.shape[0], l))
Cx_vec_xi = np.zeros((xi.shape[0]))
Cy_vec_xi = np.zeros((xi.shape[0]))
N_xi = fn.basisfunction_return_N_at_xi(N_xi, xi, knot, l, n)
for i in range(xi.shape[0]):
    Cx_xi = 0
    Cy_xi = 0
    for j in range(l):
        Cx_xi += N_xi[n][i][j] * CP[j][0]
        Cy_xi += N_xi[n][i][j] * CP[j][1]
        Cx_vec_xi[i] = Cx_xi
        Cy_vec_xi[i] = Cy_xi
ax3.scatter(Cx_vec_xi, Cy_vec_xi, c=color[2], marker="s", s=5)
ax3.set_aspect('equal', adjustable='box')
ax3.grid()
ax3.set_xlim(0, 2)
ax3.set_ylim(0, 2)


# Define control points
# # weight value
# hight_rate = math.sqrt(2) - 1.
# wv = math.cos(math.atan(hight_rate / 1.)) * math.cos(math.atan(hight_rate / 1.))

# CP_kinji_0 = np.array([[5., 0., 0., 1.],
#                        [5., 5.*hight_rate, 0., 1.*wv],
#                        [5.*hight_rate, 5., 0., 1.*wv],
#                        [0., 5., 0., 1.]])

# # CP_kinji_0 = np.array([[5., 0., 0., 1.],
# #                        [5., 5., 0., 1.*wv1],
# #                        [0., 5., 0., 1.]])

# CP_kinji = CP_kinji_0[:, :-1]
# w = CP_kinji_0[:, 3:]

# knot = np.array([0., 0., 0., .5, 1., 1., 1.])
# # knot = np.array([0., 0., 0., 1., 1., 1.])

# プログラム
# CP_kinji = np.array(
#    [[ 1.000000000000000000000e+00,   0.000000000000000000000e+00,   1.000000000000000000000e+00 ],
#     [ 1.000000000000000000000e+00,   9.752680458061155355232e-02,   9.809698831278217401319e-01 ],
#     [ 9.619397662556432582193e-01,   2.907978998723738661702e-01,   9.619397662556433692416e-01 ],
#     [ 8.858192987669302187470e-01,   4.745689648578056418771e-01,   9.619397662556433692416e-01 ],
#     [ 7.760686460529533370334e-01,   6.381449163201418084412e-01,   9.809698831278217401319e-01 ],
#     [ 7.071067811865475727373e-01,   7.071067811865475727373e-01,   1.000000000000000000000e+00 ]])
# l = 6
# n = 2
# m = l+n+1
# knot = np.array([0.000000000000000000000e+00,   0.000000000000000000000e+00,   0.000000000000000000000e+00,   2.500000000000000000000e-01,   5.000000000000000000000e-01,   7.500000000000000000000e-01,   1.000000000000000000000e+00,   1.000000000000000000000e+00,   1.000000000000000000000e+00])

# igs
# CP_kinji = np.array(
#    [[ 1.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000 ],
#     [ 1.00000000000000000000, 0.09752680458061159519, 0.98096988312782174013 ],
#     [ 0.96193976625564348026, 0.29079789987237392168, 0.96193976625564348026 ],
#     [ 0.88581929876693032977, 0.47456896485780580841, 0.96193976625564336924 ],
#     [ 0.77606864605295333703, 0.63814491632014191946, 0.98096988312782174013 ],
#     [ 0.70710678118654757274, 0.70710678118654757274, 1.00000000000000000000 ]])
# l = 6
# n = 2
# m = l+n+1
# knot = np.array([0.000000000000000000000e+00,   0.000000000000000000000e+00,   0.000000000000000000000e+00,   2.500000000000000000000e-01,   5.000000000000000000000e-01,   7.500000000000000000000e-01,   1.000000000000000000000e+00,   1.000000000000000000000e+00,   1.000000000000000000000e+00])

# 1/4円
# wv = 1. / math.sqrt(2)
# CP_kinji = np.array([[0.,1.,1.],
#                      [1.,1.,wv],
#                      [1.,0.,1.]])
# l = 3
# n = 2
# m = l+n+1
# knot = fn.def_knot(m, n)

# 1/4円 order 3
# CP_kinji = np.array([[1.0000000000000000e+00,	0.0000000000000000e+00,	1.0000000000000000e+00],
#                      [1.0000000000000000e+00,	7.8840936161031966e-02,	9.7210540773205223e-01],
#                      [9.8815430172276020e-01,	2.3691396554479624e-01,	9.2428610670128442e-01],
#                      [9.2739697290564982e-01,	4.5985470896153996e-01,	8.7646680567051660e-01],
#                      [8.1702429809208277e-01,	6.5814165201237163e-01,	8.5255715515513286e-01],
#                      [6.5814165201237163e-01,	8.1702429809208288e-01,	8.5255715515513286e-01],
#                      [4.5985470896153990e-01,	9.2739697290564971e-01,	8.7646680567051671e-01],
#                      [2.3691396554479624e-01,	9.8815430172276009e-01,	9.2428610670128442e-01],
#                      [7.8840936161031980e-02,	1.0000000000000000e+00,	9.7210540773205212e-01],
#                      [0.0000000000000000e+00,	1.0000000000000000e+00,	1.0000000000000000e+00]])

CP_kinji = np.array([[1.0000000000000000e+00,	0.0000000000000000e+00,	1.0000000000000000e+00],
                     [1.0000000000000000e+00,	6.9275926699411725e-02,	9.7210540773205223e-01],
                     [9.8528005335336311e-01,	2.1248285440492004e-01,	9.2428610670128442e-01],
                     [9.1462319908904877e-01,	4.2564697568996424e-01,	8.7646680567051660e-01],
                     [7.9254089039093678e-01,	6.2497776339900124e-01,	8.5255715515513286e-01],
                     [6.2497776339900124e-01,	7.9254089039093689e-01,	8.5255715515513286e-01],
                     [4.2564697568996418e-01,	9.1462319908904866e-01,	8.7646680567051671e-01],
                     [2.1248285440492004e-01,	9.8528005335336299e-01,	9.2428610670128442e-01],
                     [6.9275926699411738e-02,	1.0000000000000000e+00,	9.7210540773205212e-01],
                     [0.0000000000000000e+00,	1.0000000000000000e+00,	1.0000000000000000e+00]])
l = 10
n = 3
m = l+n+1
knot = fn.def_knot(m, n)


# 1/8円
# pi = math.pi
# v0 = math.cos(pi/4.)
# v1 = math.sin(pi/4.)
# CP_kinji = np.array([[1.,0.,1.],
#                      [1.,math.tan(pi/8.),math.cos(pi/8.)],
#                      [v0,v1,1.]])
# l = 3
# n = 2
# m = l+n+1
# knot = fn.def_knot(m, n)

# number_of_auto_insertion = 1
# CP_kinji, l, m, knot = fn.knot_insertion_C(CP_kinji, n, l, m, knot, number_of_auto_insertion)

# insert_knot = np.array([0.5])
# CP_kinji, l, m, knot = fn.NURBS_knot_insertion_B(CP_kinji, n, l, knot, insert_knot)

number_of_auto_insertion = 3
CP_kinji, l, m, knot = fn.NURBS_knot_insertion_C(CP_kinji, n, l, m, knot, number_of_auto_insertion)

w = CP_kinji[:,2]

# 変数宣言
delta = 500

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
    Cy_vec_2[i] = math.sin(math.acos(Cx))
    print(Cy_vec_1[i] - Cy_vec_2[i])
ax2.plot(Cx_vec_1, Cy_vec_1, c=color[3], marker="", linewidth=0.5)
ax2.plot(Cx_vec_1, Cy_vec_2, c=color[4], marker="", linewidth=1.0)

ax1.set_axisbelow(True)
ax2.set_axisbelow(True)
ax3.set_axisbelow(True)

fig.set_figheight(9)
fig.set_figwidth(12)
plt.show()
