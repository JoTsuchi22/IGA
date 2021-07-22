import numpy as np
import matplotlib.pyplot as plt
import math

# difine filename
filename = "tokuipatch"

# patch info (パッチ番号, ξ方向のコントロールポイント個数, η方向のコントロールポイント個数)
patch0_info = np.array([0, 3, 3])
patch1_info = np.array([1, 3, 3])
patch2_info = np.array([2, 3, 3])
patch3_info = np.array([3, 3, 3])
patch4_info = np.array([4, 6, 3])
patch5_info = np.array([5, 6, 6])
patch6_info = np.array([6, 3, 6])
patch7_info = np.array([7, 3, 6])
patch8_info = np.array([8, 3, 6])
patch9_info = np.array([9, 3, 3])

# パッチのコントロールポイント 座標，重み
patch0 = np.array(
   [[1.,0.,1.],
    [1.,0.,1.],
    [1.5,0.,1.],
    [1.,0.,1.],
    [1.,0.,1.],
    [1.5,0.25,1.],
    [1.,0.,1.],
    [1.,0.,1.],
    [1.5,0.5,1.]])

patch1 = np.array(
   [[1.,0.,1.],
    [1.,0.,1.],
    [1.5,0.5,1.],
    [1.,0.,1.],
    [1.,0.,1.],
    [1.25,0.5,1.],
    [1.,0.,1.],
    [1.,0.,1.],
    [1.,0.5,1.]])

patch2 = np.array(
   [[1.,0.,1.],
    [1.,0.,1.],
    [1.,0.5,1.],
    [1.,0.,1.],
    [1.,0.,1.],
    [0.75,0.5,1.],
    [1.,0.,1.],
    [1.,0.,1.],
    [0.5,0.5,1.]])

patch3 = np.array(
   [[1.,0.,1.],
    [1.,0.,1.],
    [0.5,0.5,1.],
    [1.,0.,1.],
    [1.,0.,1.],
    [0.5,0.25,1.],
    [1.,0.,1.],
    [1.,0.,1.],
    [0.5,0.,1.]])

patch4 = np.array(
   [[ 1.500000000000000000000e+00,   0.000000000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 2.562500000000000000000e+00,   0.000000000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 4.687500000000000000000e+00,   0.000000000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 6.812500000000000000000e+00,   0.000000000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 8.937500000000000000000e+00,   0.000000000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 1.000000000000000000000e+01,   0.000000000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 1.500000000000000000000e+00,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 2.562500000000000000000e+00,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 4.687500000000000000000e+00,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 6.812500000000000000000e+00,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 8.937500000000000000000e+00,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 1.000000000000000000000e+01,   5.000000000000000000000e-01,   1.000000000000000000000e+00 ],
    [ 1.500000000000000000000e+00,   1.000000000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 2.562500000000000000000e+00,   1.000000000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 4.687500000000000000000e+00,   1.000000000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 6.812500000000000000000e+00,   1.000000000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 8.937500000000000000000e+00,   1.000000000000000000000e+00,   1.000000000000000000000e+00 ],
    [ 1.000000000000000000000e+01,   1.000000000000000000000e+00,   1.000000000000000000000e+00 ]])
    
patch5 = np.array(
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

patch6 = np.array(
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

patch7 = np.array(
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

patch8 = np.array(
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

patch9 = np.array(
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
allpointnumber_overlapping = patch0_info[1] * patch0_info[2] + patch1_info[1] * patch1_info[2] + patch2_info[1] * patch2_info[2] \
                           + patch3_info[1] * patch3_info[2] + patch4_info[1] * patch4_info[2] + patch5_info[1] * patch5_info[2] \
                           + patch6_info[1] * patch6_info[2] + patch7_info[1] * patch7_info[2] + patch8_info[1] * patch8_info[2] \
                           + patch9_info[1] * patch9_info[2]

globalpoint = np.zeros((allpointnumber_overlapping))
globalpoint_x = np.zeros((allpointnumber_overlapping))
globalpoint_y = np.zeros((allpointnumber_overlapping))
globalpoint_w = np.zeros((allpointnumber_overlapping))
globalpoint_bool = np.zeros((allpointnumber_overlapping), dtype = bool)

localpoint0 = np.zeros(int(patch0_info[1] * patch0_info[2]))
localpoint1 = np.zeros(int(patch1_info[1] * patch1_info[2]))
localpoint2 = np.zeros(int(patch2_info[1] * patch2_info[2]))
localpoint3 = np.zeros(int(patch3_info[1] * patch3_info[2]))
localpoint4 = np.zeros(int(patch4_info[1] * patch4_info[2]))
localpoint5 = np.zeros(int(patch5_info[1] * patch5_info[2]))
localpoint6 = np.zeros(int(patch6_info[1] * patch6_info[2]))
localpoint7 = np.zeros(int(patch7_info[1] * patch7_info[2]))
localpoint8 = np.zeros(int(patch8_info[1] * patch8_info[2]))
localpoint9 = np.zeros(int(patch9_info[1] * patch9_info[2]))

A0_i = np.zeros((2, patch0_info[2]))
A0_j = np.zeros((2, patch0_info[1]))
A1_i = np.zeros((2, patch1_info[2]))
A1_j = np.zeros((2, patch1_info[1]))
A2_i = np.zeros((2, patch2_info[2]))
A2_j = np.zeros((2, patch2_info[1]))
A3_i = np.zeros((2, patch3_info[2]))
A3_j = np.zeros((2, patch3_info[1]))
A4_i = np.zeros((2, patch4_info[2]))
A4_j = np.zeros((2, patch4_info[1]))
A5_i = np.zeros((2, patch5_info[2]))
A5_j = np.zeros((2, patch5_info[1]))
A6_i = np.zeros((2, patch6_info[2]))
A6_j = np.zeros((2, patch6_info[1]))
A7_i = np.zeros((2, patch7_info[2]))
A7_j = np.zeros((2, patch7_info[1]))
A8_i = np.zeros((2, patch8_info[2]))
A8_j = np.zeros((2, patch8_info[1]))
A9_i = np.zeros((2, patch9_info[2]))
A9_j = np.zeros((2, patch9_info[1]))

# globalpatch = patch0
length = int(globalpoint[globalpoint_bool].shape[0])
a_g = length
a_l_g = 0
a_l_l = 0
patch0_x = np.reshape(patch0[:,0], [patch0_info[1], patch0_info[2]]).T
patch0_y = np.reshape(patch0[:,1], [patch0_info[1], patch0_info[2]]).T
patch0_w = np.reshape(patch0[:,2], [patch0_info[1], patch0_info[2]]).T
for j in range(patch0_info[2]):
    for i in range(patch0_info[1]):
        localpoint0[a_l_g] = length + a_l_l
        globalpoint[a_g] = a_g
        globalpoint_bool[a_g] = True
        globalpoint_x[a_g] = patch0_x[i][j]
        globalpoint_y[a_g] = patch0_y[i][j]
        globalpoint_w[a_g] = patch0_w[i][j]
        a_g = a_g + 1
        a_l_g = a_l_g + 1
        a_l_l = a_l_l + 1
a_l_g = 0
for j in range(patch0_info[2]):
    for i in range(patch0_info[1]):
        if i == 0:
            A0_i[0,j] = localpoint0[a_l_g]
        if i == patch0_info[1] - 1:
            A0_i[1,j] = localpoint0[a_l_g]
        if j == 0:
            A0_j[0,i] = localpoint0[a_l_g]
        if j == patch0_info[2] - 1:
            A0_j[1,i] = localpoint0[a_l_g]
        a_l_g = a_l_g + 1

# globalpatch += patch1
length = int(globalpoint[globalpoint_bool].shape[0])
a_g = length
a_l_g = 0
a_l_l = 0
patch1_x = np.reshape(patch1[:,0], [patch1_info[1], patch1_info[2]]).T
patch1_y = np.reshape(patch1[:,1], [patch1_info[1], patch1_info[2]]).T
patch1_w = np.reshape(patch1[:,2], [patch1_info[1], patch1_info[2]]).T
for j in range(patch1_info[2]):
    for i in range(patch1_info[1]):
        if j == 0:
            localpoint1[a_l_g] = A0_j[1,i]
            a_l_g = a_l_g + 1
        else:
            localpoint1[a_l_g] = length + a_l_l
            globalpoint[a_g] = a_g
            globalpoint_bool[a_g] = True
            globalpoint_x[a_g] = patch1_x[i][j]
            globalpoint_y[a_g] = patch1_y[i][j]
            globalpoint_w[a_g] = patch1_w[i][j]
            a_g = a_g + 1
            a_l_g = a_l_g + 1
            a_l_l = a_l_l + 1
a_l_g = 0
for j in range(patch1_info[2]):
    for i in range(patch1_info[1]):
        if i == 0:
            A1_i[0,j] = localpoint1[a_l_g]
        if i == patch1_info[1] - 1:
            A1_i[1,j] = localpoint1[a_l_g]
        if j == 0:
            A1_j[0,i] = localpoint1[a_l_g]
        if j == patch1_info[2] - 1:
            A1_j[1,i] = localpoint1[a_l_g]
        a_l_g = a_l_g + 1

# globalpatch += patch2
length = int(globalpoint[globalpoint_bool].shape[0])
a_g = length
a_l_g = 0
a_l_l = 0
patch2_x = np.reshape(patch2[:,0], [patch2_info[1], patch2_info[2]]).T
patch2_y = np.reshape(patch2[:,1], [patch2_info[1], patch2_info[2]]).T
patch2_w = np.reshape(patch2[:,2], [patch2_info[1], patch2_info[2]]).T
for j in range(patch2_info[2]):
    for i in range(patch2_info[1]):
        if j == 0:
            localpoint2[a_l_g] = A1_j[1,i]
            a_l_g = a_l_g + 1
        else:
            localpoint2[a_l_g] = length + a_l_l
            globalpoint[a_g] = a_g
            globalpoint_bool[a_g] = True
            globalpoint_x[a_g] = patch2_x[i][j]
            globalpoint_y[a_g] = patch2_y[i][j]
            globalpoint_w[a_g] = patch2_w[i][j]
            a_g = a_g + 1
            a_l_g = a_l_g + 1
            a_l_l = a_l_l + 1
a_l_g = 0
for j in range(patch2_info[2]):
    for i in range(patch2_info[1]):
        if i == 0:
            A2_i[0,j] = localpoint2[a_l_g]
        if i == patch2_info[1] - 1:
            A2_i[1,j] = localpoint2[a_l_g]
        if j == 0:
            A2_j[0,i] = localpoint2[a_l_g]
        if j == patch2_info[2] - 1:
            A2_j[1,i] = localpoint2[a_l_g]
        a_l_g = a_l_g + 1

# globalpatch += patch3
length = int(globalpoint[globalpoint_bool].shape[0])
a_g = length
a_l_g = 0
a_l_l = 0
patch3_x = np.reshape(patch3[:,0], [patch3_info[1], patch3_info[2]]).T
patch3_y = np.reshape(patch3[:,1], [patch3_info[1], patch3_info[2]]).T
patch3_w = np.reshape(patch3[:,2], [patch3_info[1], patch3_info[2]]).T
for j in range(patch3_info[2]):
    for i in range(patch3_info[1]):
        if j == 0:
            localpoint3[a_l_g] = A2_j[1,i]
            a_l_g = a_l_g + 1
        else:
            localpoint3[a_l_g] = length + a_l_l
            globalpoint[a_g] = a_g
            globalpoint_bool[a_g] = True
            globalpoint_x[a_g] = patch3_x[i][j]
            globalpoint_y[a_g] = patch3_y[i][j]
            globalpoint_w[a_g] = patch3_w[i][j]
            a_g = a_g + 1
            a_l_g = a_l_g + 1
            a_l_l = a_l_l + 1
a_l_g = 0
for j in range(patch3_info[2]):
    for i in range(patch3_info[1]):
        if i == 0:
            A3_i[0,j] = localpoint3[a_l_g]
        if i == patch3_info[1] - 1:
            A3_i[1,j] = localpoint3[a_l_g]
        if j == 0:
            A3_j[0,i] = localpoint3[a_l_g]
        if j == patch3_info[2] - 1:
            A3_j[1,i] = localpoint3[a_l_g]
        a_l_g = a_l_g + 1

# # globalpatch += patch4
# length = int(globalpoint[globalpoint_bool].shape[0])
# a_g = length
# a_l_g = 0
# a_l_l = 0
# patch4_x = np.reshape(patch4[:,0], [patch4_info[1], patch4_info[2]]).T
# patch4_y = np.reshape(patch4[:,1], [patch4_info[1], patch4_info[2]]).T
# patch4_w = np.reshape(patch4[:,2], [patch4_info[1], patch4_info[2]]).T
# for j in range(patch4_info[2]):
#     for i in range(patch4_info[1]):
#         if i == 0:
#             localpoint4[a_l_g] = A0_i[1,j]
#             a_l_g = a_l_g + 1
#         else:
#             localpoint4[a_l_g] = length + a_l_l
#             globalpoint[a_g] = a_g
#             globalpoint_bool[a_g] = True
#             globalpoint_x[a_g] = patch4_x[i][j]
#             globalpoint_y[a_g] = patch4_y[i][j]
#             globalpoint_w[a_g] = patch4_w[i][j]
#             a_g = a_g + 1
#             a_l_g = a_l_g + 1
#             a_l_l = a_l_l + 1
# a_l_g = 0
# for j in range(patch4_info[2]):
#     for i in range(patch4_info[1]):
#         if i == 0:
#             A4_i[0,j] = localpoint4[a_l_g]
#         if i == patch4_info[1] - 1:
#             A4_i[1,j] = localpoint4[a_l_g]
#         if j == 0:
#             A4_j[0,i] = localpoint4[a_l_g]
#         if j == patch4_info[2] - 1:
#             A4_j[1,i] = localpoint4[a_l_g]
#         a_l_g = a_l_g + 1

# # globalpatch += patch5
# length = int(globalpoint[globalpoint_bool].shape[0])
# a_g = length
# a_l_g = 0
# a_l_l = 0
# patch5_x = np.reshape(patch5[:,0], [patch5_info[1], patch5_info[2]]).T
# patch5_y = np.reshape(patch5[:,1], [patch5_info[1], patch5_info[2]]).T
# patch5_w = np.reshape(patch5[:,2], [patch5_info[1], patch5_info[2]]).T
# for j in range(patch5_info[2]):
#     for i in range(patch5_info[1]):
#         if j == 0:
#             localpoint5[a_l_g] = A4_j[1,i]
#             a_l_g = a_l_g + 1
#         else:
#             localpoint5[a_l_g] = length + a_l_l
#             globalpoint[a_g] = a_g
#             globalpoint_bool[a_g] = True
#             globalpoint_x[a_g] = patch5_x[i][j]
#             globalpoint_y[a_g] = patch5_y[i][j]
#             globalpoint_w[a_g] = patch5_w[i][j]
#             a_g = a_g + 1
#             a_l_g = a_l_g + 1
#             a_l_l = a_l_l + 1
# a_l_g = 0
# for j in range(patch5_info[2]):
#     for i in range(patch5_info[1]):
#         if i == 0:
#             A5_i[0,j] = localpoint5[a_l_g]
#         if i == patch5_info[1] - 1:
#             A5_i[1,j] = localpoint5[a_l_g]
#         if j == 0:
#             A5_j[0,i] = localpoint5[a_l_g]
#         if j == patch5_info[2] - 1:
#             A5_j[1,i] = localpoint5[a_l_g]
#         a_l_g = a_l_g + 1

# # globalpatch += patch6
# length = int(globalpoint[globalpoint_bool].shape[0])
# a_g = length
# a_l_g = 0
# a_l_l = 0
# patch6_x = np.reshape(patch6[:,0], [patch6_info[1], patch6_info[2]]).T
# patch6_y = np.reshape(patch6[:,1], [patch6_info[1], patch6_info[2]]).T
# patch6_w = np.reshape(patch6[:,2], [patch6_info[1], patch6_info[2]]).T
# for j in range(patch6_info[2]):
#     for i in range(patch6_info[1]):
#         if i == patch6_info[1]-1:
#             localpoint6[a_l_g] = A5_i[0,j]
#             a_l_g = a_l_g + 1
#         else:
#             localpoint6[a_l_g] = length + a_l_l
#             globalpoint[a_g] = a_g
#             globalpoint_bool[a_g] = True
#             globalpoint_x[a_g] = patch6_x[i][j]
#             globalpoint_y[a_g] = patch6_y[i][j]
#             globalpoint_w[a_g] = patch6_w[i][j]
#             a_g = a_g + 1
#             a_l_g = a_l_g + 1
#             a_l_l = a_l_l + 1
# a_l_g = 0
# for j in range(patch6_info[2]):
#     for i in range(patch6_info[1]):
#         if i == 0:
#             A6_i[0,j] = localpoint6[a_l_g]
#         if i == patch6_info[1] - 1:
#             A6_i[1,j] = localpoint6[a_l_g]
#         if j == 0:
#             A6_j[0,i] = localpoint6[a_l_g]
#         if j == patch6_info[2] - 1:
#             A6_j[1,i] = localpoint6[a_l_g]
#         a_l_g = a_l_g + 1

# # globalpatch += patch7
# length = int(globalpoint[globalpoint_bool].shape[0])
# a_g = length
# a_l_g = 0
# a_l_l = 0
# patch7_x = np.reshape(patch7[:,0], [patch7_info[1], patch7_info[2]]).T
# patch7_y = np.reshape(patch7[:,1], [patch7_info[1], patch7_info[2]]).T
# patch7_w = np.reshape(patch7[:,2], [patch7_info[1], patch7_info[2]]).T
# for j in range(patch7_info[2]):
#     for i in range(patch7_info[1]):
#         if i == patch7_info[1]-1:
#             localpoint7[a_l_g] = A6_i[0,j]
#             a_l_g = a_l_g + 1
#         else:
#             localpoint7[a_l_g] = length + a_l_l
#             globalpoint[a_g] = a_g
#             globalpoint_bool[a_g] = True
#             globalpoint_x[a_g] = patch7_x[i][j]
#             globalpoint_y[a_g] = patch7_y[i][j]
#             globalpoint_w[a_g] = patch7_w[i][j]
#             a_g = a_g + 1
#             a_l_g = a_l_g + 1
#             a_l_l = a_l_l + 1
# a_l_g = 0
# for j in range(patch7_info[2]):
#     for i in range(patch7_info[1]):
#         if i == 0:
#             A7_i[0,j] = localpoint7[a_l_g]
#         if i == patch7_info[1] - 1:
#             A7_i[1,j] = localpoint7[a_l_g]
#         if j == 0:
#             A7_j[0,i] = localpoint7[a_l_g]
#         if j == patch7_info[2] - 1:
#             A7_j[1,i] = localpoint7[a_l_g]
#         a_l_g = a_l_g + 1

# # globalpatch += patch8
# length = int(globalpoint[globalpoint_bool].shape[0])
# a_g = length
# a_l_g = 0
# a_l_l = 0
# patch8_x = np.reshape(patch8[:,0], [patch8_info[1], patch8_info[2]]).T
# patch8_y = np.reshape(patch8[:,1], [patch8_info[1], patch8_info[2]]).T
# patch8_w = np.reshape(patch8[:,2], [patch8_info[1], patch8_info[2]]).T
# for j in range(patch8_info[2]):
#     for i in range(patch8_info[1]):
#         if i == patch8_info[1]-1:
#             localpoint8[a_l_g] = A7_i[0,j]
#             a_l_g = a_l_g + 1
#         else:
#             localpoint8[a_l_g] = length + a_l_l
#             globalpoint[a_g] = a_g
#             globalpoint_bool[a_g] = True
#             globalpoint_x[a_g] = patch8_x[i][j]
#             globalpoint_y[a_g] = patch8_y[i][j]
#             globalpoint_w[a_g] = patch8_w[i][j]
#             a_g = a_g + 1
#             a_l_g = a_l_g + 1
#             a_l_l = a_l_l + 1
# a_l_g = 0
# for j in range(patch8_info[2]):
#     for i in range(patch8_info[1]):
#         if i == 0:
#             A8_i[0,j] = localpoint8[a_l_g]
#         if i == patch8_info[1] - 1:
#             A8_i[1,j] = localpoint8[a_l_g]
#         if j == 0:
#             A8_j[0,i] = localpoint8[a_l_g]
#         if j == patch8_info[2] - 1:
#             A8_j[1,i] = localpoint8[a_l_g]
#         a_l_g = a_l_g + 1




# # globalpatch += patch9
# length = int(globalpoint[globalpoint_bool].shape[0])
# a_g = length
# a_l_g = 0
# a_l_l = 0
# patch9_x = np.reshape(patch9[:,0], [patch9_info[1], patch9_info[2]]).T
# patch9_y = np.reshape(patch9[:,1], [patch9_info[1], patch9_info[2]]).T
# patch9_w = np.reshape(patch9[:,2], [patch9_info[1], patch9_info[2]]).T
# for j in range(patch9_info[2]):
#     for i in range(patch9_info[1]):
#         if i == patch9_info[1] - 1 and j == patch9_info[2] - 1:
#             localpoint9[a_l_g] = A3_j[1,patch3_info[2]-1-j]
#             a_l_g = a_l_g + 1
#         if i == patch9_info[1] - 1 and j != patch9_info[2] - 1:
#             localpoint9[a_l_g] = A3_j[1,patch3_info[2]-1-j]
#             a_l_g = a_l_g + 1
#         if i != patch9_info[1] - 1 and j == patch9_info[2] - 1:
#             localpoint9[a_l_g] = A8_j[0,i]
#             a_l_g = a_l_g + 1
#         if i != patch9_info[1] - 1 and j != patch9_info[2] - 1:
#             localpoint9[a_l_g] = length + a_l_l
#             globalpoint[a_g] = a_g
#             globalpoint_bool[a_g] = True
#             globalpoint_x[a_g] = patch9_x[i][j]
#             globalpoint_y[a_g] = patch9_y[i][j]
#             globalpoint_w[a_g] = patch9_w[i][j]
#             a_g = a_g + 1
#             a_l_g = a_l_g + 1
#             a_l_l = a_l_l + 1
# a_l_g = 0
# for j in range(patch9_info[2]):
#     for i in range(patch9_info[1]):
#         if i == 0:
#             A9_i[0,j] = localpoint9[a_l_g]
#         if i == patch1_info[1] - 1:
#             A9_i[1,j] = localpoint9[a_l_g]
#         if j == 0:
#             A9_j[0,i] = localpoint9[a_l_g]
#         if j == patch1_info[2] - 1:
#             A9_j[1,i] = localpoint9[a_l_g]
#         a_l_g = a_l_g + 1


# テキストデータ書き込み
def write_date_header(filename):
    filename_txt = filename + ".txt"
    f = open(filename_txt, 'w')
    f.write('patch connectivity')
    f.write('\n')

def write_date_localpoint(filename, localpoint):
    filename_txt = filename + ".txt"
    f = open(filename_txt, 'a')
    for i in range(localpoint.shape[0]):
        if i == 0:
            f.write(str(int(localpoint[i])))
        if i != 0 and i != localpoint.shape[0] - 1:
            f.write('   ')
            f.write(str(int(localpoint[i])))
        if i == localpoint.shape[0] - 1:
            f.write('   ')
            f.write(str(int(localpoint[i])))
            f.write('\n')

def write_date_globalpoint(filename, globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w):
    filename_txt = filename + ".txt"
    f = open(filename_txt, 'a')
    f.write('\n')
    A = np.zeros((globalpoint[globalpoint_bool].shape[0],4))
    for i in range(globalpoint[globalpoint_bool].shape[0]):
        A[i][0] = globalpoint[globalpoint_bool][i]
        A[i][1] = globalpoint_x[globalpoint_bool][i]
        A[i][2] = globalpoint_y[globalpoint_bool][i]
        A[i][3] = globalpoint_w[globalpoint_bool][i]
        f.write(str(int(A[i][0])))
        f.write('   ')
        f.write('   ')
        f.write(str('{:.21e}'.format(A[i][1])))
        f.write('   ')
        f.write('   ')
        f.write(str('{:.21e}'.format(A[i][2])))
        f.write('   ')
        f.write('   ')
        f.write(str('{:.21e}'.format(A[i][3])))
        f.write('\n')
    f.close()


write_date_header(filename)
write_date_localpoint(filename, localpoint0)
write_date_localpoint(filename, localpoint1)
write_date_localpoint(filename, localpoint2)
write_date_localpoint(filename, localpoint3)
write_date_localpoint(filename, localpoint4)
write_date_globalpoint(filename, globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w)

def write_BC(filename, A, axis):
    filename_txt = filename + ".txt"
    f = open(filename_txt, 'a')
    f.write('\n')
    for i in range(6):
        f.write(str(int(A[axis,i])))
        f.write('\n')
    f.write('\n')

# 境界条件txt出力
filename = filename + '_BC'
filename_txt = filename + ".txt"
f = open(filename_txt, 'w')

write_BC(filename, A0_j, 0)
write_BC(filename, A1_j, 0)
write_BC(filename, A3_i, 0)
write_BC(filename, A4_j, 1)
write_BC(filename, A2_j, 1)
write_BC(filename, A3_j, 1)


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