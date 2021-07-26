import numpy as np
import matplotlib.pyplot as plt
import math
import make_inputdate_function as mif

# difine filename
filename = "tokuipatch"

# patch info (ξ方向のコントロールポイント個数, η方向のコントロールポイント個数)
patch_info = np.array([[3, 3],
                       [3, 3],
                       [3, 3],
                       [3, 3],
                       [6, 3],
                       [6, 6],
                       [3, 6],
                       [3, 6],
                       [3, 6],
                       [3, 3]])

xi_max = 0
eta_max = 0
for i in range(patch_info.shape[0]):
    if xi_max < patch_info[i][0]:
        xi_max = patch_info[i][0]
    if eta_max < patch_info[i][1]:
        eta_max = patch_info[i][1]

# パッチのコントロールポイント 座標，重み
patch = np.zeros((patch_info.shape[0], xi_max * eta_max, 3))
patch_bool = np.zeros((patch_info.shape[0], xi_max * eta_max, 3), dtype=bool)
for i in range(patch_info.shape[0]):
    patch_bool[i,:patch_info[i][0]*patch_info[i][1],:] = True

# 各パッチ入力
num = 0
patch[num,:patch_info[num][0]*patch_info[num][1],:] = np.array(
   [[1.,0.,1.],
    [1.,0.,1.],
    [1.5,0.,1.],
    [1.,0.,1.],
    [1.,0.,1.],
    [1.5,0.25,1.],
    [1.,0.,1.],
    [1.,0.,1.],
    [1.5,0.5,1.]])

num = 1
patch[num,:patch_info[num][0]*patch_info[num][1],:] = np.array(
   [[1.,0.,1.],
    [1.,0.,1.],
    [1.5,0.5,1.],
    [1.,0.,1.],
    [1.,0.,1.],
    [1.25,0.5,1.],
    [1.,0.,1.],
    [1.,0.,1.],
    [1.,0.5,1.]])

num = 2
patch[num,:patch_info[num][0]*patch_info[num][1],:] = np.array(
   [[1.,0.,1.],
    [1.,0.,1.],
    [1.,0.5,1.],
    [1.,0.,1.],
    [1.,0.,1.],
    [0.75,0.5,1.],
    [1.,0.,1.],
    [1.,0.,1.],
    [0.5,0.5,1.]])

num = 3
patch[num,:patch_info[num][0]*patch_info[num][1],:] = np.array(
   [[1.,0.,1.],
    [1.,0.,1.],
    [0.5,0.5,1.],
    [1.,0.,1.],
    [1.,0.,1.],
    [0.5,0.25,1.],
    [1.,0.,1.],
    [1.,0.,1.],
    [0.5,0.,1.]])

num = 4
patch[num,:patch_info[num][0]*patch_info[num][1],:] = np.array(
   [[ 1.500000000000000000000e+00,   0.000000000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 2.562500000000000000000e+00,   0.000000000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 4.687500000000000000000e+00,   0.000000000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 6.812500000000000000000e+00,   0.000000000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 8.937500000000000000000e+00,   0.000000000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 1.000000000000000000000e+01,   0.000000000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 1.500000000000000000000e+00,   2.500000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 2.562500000000000000000e+00,   2.500000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 4.687500000000000000000e+00,   2.500000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 6.812500000000000000000e+00,   2.500000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 8.937500000000000000000e+00,   2.500000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 1.000000000000000000000e+01,   2.500000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 1.500000000000000000000e+00,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 2.562500000000000000000e+00,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 4.687500000000000000000e+00,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 6.812500000000000000000e+00,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 8.937500000000000000000e+00,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 1.000000000000000000000e+01,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ]])

num = 5
patch[num,:patch_info[num][0]*patch_info[num][1],:] = np.array(
   [[ 1.500000000000000000000e+00,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 2.562500000000000000000e+00,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 4.687500000000000000000e+00,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 6.812500000000000000000e+00,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 8.937500000000000000000e+00,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 1.000000000000000000000e+01,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 1.500000000000000000000e+00,   1.687500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 2.562500000000000000000e+00,   1.687500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 4.687500000000000000000e+00,   1.687500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 6.812500000000000000000e+00,   1.687500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 8.937500000000000000000e+00,   1.687500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 1.000000000000000000000e+01,   1.687500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 1.500000000000000000000e+00,   4.062500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 2.562500000000000000000e+00,   4.062500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 4.687500000000000000000e+00,   4.062500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 6.812500000000000000000e+00,   4.062500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 8.937500000000000000000e+00,   4.062500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 1.000000000000000000000e+01,   4.062500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 1.500000000000000000000e+00,   6.437500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 2.562500000000000000000e+00,   6.437500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 4.687500000000000000000e+00,   6.437500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 6.812500000000000000000e+00,   6.437500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 8.937500000000000000000e+00,   6.437500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 1.000000000000000000000e+01,   6.437500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 1.500000000000000000000e+00,   8.812500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 2.562500000000000000000e+00,   8.812500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 4.687500000000000000000e+00,   8.812500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 6.812500000000000000000e+00,   8.812500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 8.937500000000000000000e+00,   8.812500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 1.000000000000000000000e+01,   8.812500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 1.500000000000000000000e+00,   1.000000000000000000000e+01,   1.000000000000000000000e+00 ],
    [ 2.562500000000000000000e+00,   1.000000000000000000000e+01,   1.000000000000000000000e+00 ],
    [ 4.687500000000000000000e+00,   1.000000000000000000000e+01,   1.000000000000000000000e+00 ],
    [ 6.812500000000000000000e+00,   1.000000000000000000000e+01,   1.000000000000000000000e+00 ],
    [ 8.937500000000000000000e+00,   1.000000000000000000000e+01,   1.000000000000000000000e+00 ],
    [ 1.000000000000000000000e+01,   1.000000000000000000000e+01,   1.000000000000000000000e+00 ]])

num = 6
patch[num,:patch_info[num][0]*patch_info[num][1],:] = np.array(
   [[ 1.000000000000000000000e+00,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 1.250000000000000000000e+00,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 1.500000000000000000000e+00,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 1.000000000000000000000e+00,   1.687500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 1.250000000000000000000e+00,   1.687500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 1.500000000000000000000e+00,   1.687500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 1.000000000000000000000e+00,   4.062500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 1.250000000000000000000e+00,   4.062500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 1.500000000000000000000e+00,   4.062500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 1.000000000000000000000e+00,   6.437500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 1.250000000000000000000e+00,   6.437500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 1.500000000000000000000e+00,   6.437500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 1.000000000000000000000e+00,   8.812500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 1.250000000000000000000e+00,   8.812500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 1.500000000000000000000e+00,   8.812500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 1.000000000000000000000e+00,   1.000000000000000000000e+01,   1.000000000000000000000e+00 ],
    [ 1.250000000000000000000e+00,   1.000000000000000000000e+01,   1.000000000000000000000e+00 ],
    [ 1.500000000000000000000e+00,   1.000000000000000000000e+01,   1.000000000000000000000e+00 ]])

num = 7
patch[num,:patch_info[num][0]*patch_info[num][1],:] = np.array(
   [[ 5.000000000000000000000e-01,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 7.500000000000000000000e-01,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 1.000000000000000000000e+00,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 5.000000000000000000000e-01,   1.687500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 7.500000000000000000000e-01,   1.687500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 1.000000000000000000000e+00,   1.687500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 5.000000000000000000000e-01,   4.062500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 7.500000000000000000000e-01,   4.062500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 1.000000000000000000000e+00,   4.062500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 5.000000000000000000000e-01,   6.437500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 7.500000000000000000000e-01,   6.437500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 1.000000000000000000000e+00,   6.437500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 5.000000000000000000000e-01,   8.812500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 7.500000000000000000000e-01,   8.812500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 1.000000000000000000000e+00,   8.812500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 5.000000000000000000000e-01,   1.000000000000000000000e+01,   1.000000000000000000000e+00 ],
    [ 7.500000000000000000000e-01,   1.000000000000000000000e+01,   1.000000000000000000000e+00 ],
    [ 1.000000000000000000000e+00,   1.000000000000000000000e+01,   1.000000000000000000000e+00 ]])

num = 8
patch[num,:patch_info[num][0]*patch_info[num][1],:] = np.array(
   [[ 0.000000000000000000000e+00,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 2.500000000000000000000e-01,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 5.000000000000000000000e-01,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 0.000000000000000000000e+00,   1.687500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 2.500000000000000000000e-01,   1.687500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 5.000000000000000000000e-01,   1.687500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 0.000000000000000000000e+00,   4.062500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 2.500000000000000000000e-01,   4.062500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 5.000000000000000000000e-01,   4.062500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 0.000000000000000000000e+00,   6.437500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 2.500000000000000000000e-01,   6.437500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 5.000000000000000000000e-01,   6.437500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 0.000000000000000000000e+00,   8.812500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 2.500000000000000000000e-01,   8.812500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 5.000000000000000000000e-01,   8.812500000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 0.000000000000000000000e+00,   1.000000000000000000000e+01,   1.000000000000000000000e+00 ],
    [ 2.500000000000000000000e-01,   1.000000000000000000000e+01,   1.000000000000000000000e+00 ],
    [ 5.000000000000000000000e-01,   1.000000000000000000000e+01,   1.000000000000000000000e+00 ]])

num = 9
patch[num,:patch_info[num][0]*patch_info[num][1],:] = np.array(
   [[ 0.000000000000000000000e+00,   0.000000000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 2.500000000000000000000e-01,   0.000000000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 5.000000000000000000000e-01,   0.000000000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 0.000000000000000000000e+00,   2.500000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 2.500000000000000000000e-01,   2.500000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 5.000000000000000000000e-01,   2.500000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 0.000000000000000000000e+00,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 2.500000000000000000000e-01,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 5.000000000000000000000e-01,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ]])


# 重複ありの総コントロールポイント数
allpointnumber_overlapping = 0
for i in range(patch_info.shape[0]):
    allpointnumber_overlapping += patch_info[i][0] * patch_info[i][1]

globalpoint = np.zeros((allpointnumber_overlapping))
globalpoint_x = np.zeros((allpointnumber_overlapping))
globalpoint_y = np.zeros((allpointnumber_overlapping))
globalpoint_w = np.zeros((allpointnumber_overlapping))
globalpoint_bool = np.zeros((allpointnumber_overlapping), dtype = bool)

localpoint = np.zeros((patch_info.shape[0], xi_max * eta_max))
localpoint_bool = np.zeros((patch_info.shape[0], xi_max * eta_max), dtype=bool)
for i in range(patch_info.shape[0]):
    localpoint_bool[i,:patch_info[i][0]*patch_info[i][1]] = True

xi_and_eta_max = xi_max
if xi_max < eta_max:
    xi_and_eta_max = eta_max
A = np.zeros((patch_info.shape[0], 2, 2, xi_and_eta_max))
A_bool = np.zeros((patch_info.shape[0], 2, 2, xi_and_eta_max), dtype=bool)

# connect patch
xi = 0
eta = 1
positive = 0
negative = 1

patch_number = 0
globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool = \
    mif.connect_patch_arg_0boundary(patch_number, patch_info, patch, patch_bool, globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool)

patch_number = 1
connect_patch_number = 0
connection_vec = np.array([patch_number, eta, 0, connect_patch_number, eta, 1, positive])
globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool = \
    mif.connect_patch_arg_1boundary(patch_number, connection_vec, patch_info, patch, patch_bool, globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool)

patch_number = 2
connect_patch_number = 1
connection_vec = np.array([patch_number, eta, 0, connect_patch_number, eta, 1, positive])
globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool = \
    mif.connect_patch_arg_1boundary(patch_number, connection_vec, patch_info, patch, patch_bool, globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool)

patch_number = 3
connect_patch_number = 2
connection_vec = np.array([patch_number, eta, 0, connect_patch_number, eta, 1, positive])
globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool = \
    mif.connect_patch_arg_1boundary(patch_number, connection_vec, patch_info, patch, patch_bool, globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool)

patch_number = 4
connect_patch_number = 0
connection_vec = np.array([patch_number, xi, 0, connect_patch_number, xi, 1, positive])
globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool = \
    mif.connect_patch_arg_1boundary(patch_number, connection_vec, patch_info, patch, patch_bool, globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool)

patch_number = 5
connect_patch_number = 4
connection_vec = np.array([patch_number, eta, 0, connect_patch_number, eta, 1, positive])
globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool = \
    mif.connect_patch_arg_1boundary(patch_number, connection_vec, patch_info, patch, patch_bool, globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool)

patch_number = 6
connect_patch_number_1 = 1
connect_patch_number_2 = 5
connection_vec = np.array([[patch_number, eta, 0, connect_patch_number_1, xi, 1, negative],
                           [patch_number, xi, 1, connect_patch_number_2, xi, 0, positive]])
globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool = \
    mif.connect_patch_arg_2boundary(patch_number, connection_vec, patch_info, patch, patch_bool, globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool)

patch_number = 7
connect_patch_number_1 = 2
connect_patch_number_2 = 6
connection_vec = np.array([[patch_number, eta, 0, connect_patch_number_1, xi, 1, negative],
                           [patch_number, xi, 1, connect_patch_number_2, xi, 0, positive]])
globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool = \
    mif.connect_patch_arg_2boundary(patch_number, connection_vec, patch_info, patch, patch_bool, globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool)

patch_number = 8
connect_patch_number = 7
connection_vec = np.array([patch_number, xi, 1, connect_patch_number, xi, 0, positive])
globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool = \
    mif.connect_patch_arg_1boundary(patch_number, connection_vec, patch_info, patch, patch_bool, globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool)

patch_number = 9
connect_patch_number_1 = 3
connect_patch_number_2 = 8
connection_vec = np.array([[patch_number, xi, 1, connect_patch_number_1, xi, 1, negative],
                           [patch_number, eta, 1, connect_patch_number_2, eta, 0, positive]])
globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool = \
    mif.connect_patch_arg_2boundary(patch_number, connection_vec, patch_info, patch, patch_bool, globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool)


# 座標，パッチコネクティビティー，境界の辺をtxt出力
# [patch_number, xi_or_eta, 0_or_1(int)] (auto marge)
boundary_array_0 = np.array([[0, eta, 0],
                             [4, eta, 0],
                             [0,  xi, 0],
                             [1,  xi, 0],
                             [2,  xi, 0],
                             [3,  xi, 0]])

boundary_array_1 = np.array([[8,  xi, 0],
                             [9,  xi, 0]])

boundary_number = 0
mif.write_date_header(filename)
mif.write_date_localpoint(filename, patch_info, localpoint, localpoint_bool)
mif.write_date_globalpoint(filename, globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w)
boundary_number = mif.write_boundary(filename, A, A_bool, boundary_array_0, boundary_number)
boundary_number = mif.write_boundary(filename, A, A_bool, boundary_array_1, boundary_number)


# 描写
color = np.array(["r", "g", "b", "c", "m", "y", "k"])
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)

A = np.zeros((globalpoint[globalpoint_bool].shape[0],3))
for i in range(globalpoint[globalpoint_bool].shape[0]):
    A[i][0] = globalpoint[globalpoint_bool][i]
    A[i][1] = globalpoint_x[globalpoint_bool][i]
    A[i][2] = globalpoint_y[globalpoint_bool][i]
ax1.plot(A[:,1], A[:,2], c=color[0], marker="", linewidth=0.7)
for i in range(globalpoint[globalpoint_bool].shape[0]):
    ax1.text(A[i,1], A[i,2], str(int(A[i,0])), c=color[3], fontsize=6)


ax1.set_aspect('equal', adjustable='box')
ax1.set_axisbelow(True)
fig.set_figheight(9)
fig.set_figwidth(12)
ax1.grid()
ax1.set_xlim(-1, 11)
ax1.set_ylim(-1, 11)
plt.show()