import numpy as np
import matplotlib.pyplot as plt
import math
import make_inputdate_function as mif

# difine filename
filename = "connectivity_arc_loc3_11x11"

# patch info (ξ方向のコントロールポイント個数, η方向のコントロールポイント個数)
patch_info = np.array([11, 11])

# patch info 1行2列の時の例外処理
if np.array(patch_info.shape).shape[0] == 1:
    patch_info = np.array([[patch_info[0], patch_info[1]]])

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
   [[1.0000000000000000e+00,	0.0000000000000000e+00,	1.0000000000000000e+00],
    [1.0416666666666667e+00,	0.0000000000000000e+00,	1.0000000000000000e+00],
    [1.1249999999999998e+00,	0.0000000000000000e+00,	1.0000000000000000e+00],
    [1.2500000000000000e+00,	0.0000000000000000e+00,	1.0000000000000000e+00],
    [1.3750000000000000e+00,	0.0000000000000000e+00,	1.0000000000000000e+00],
    [1.5000000000000000e+00,	0.0000000000000000e+00,	1.0000000000000000e+00],
    [1.6250000000000000e+00,	0.0000000000000000e+00,	1.0000000000000000e+00],
    [1.7500000000000000e+00,	0.0000000000000000e+00,	1.0000000000000000e+00],
    [1.8750000000000000e+00,	0.0000000000000000e+00,	1.0000000000000000e+00],
    [1.9583333333333333e+00,	0.0000000000000000e+00,	1.0000000000000000e+00],
    [2.0000000000000000e+00,	0.0000000000000000e+00,	1.0000000000000000e+00],
    [1.0000000000000000e+00,	6.0399789153958709e-02,	9.7559223176554566e-01],
    [1.0416666666666667e+00,	6.2916447035373649e-02,	9.7559223176554566e-01],
    [1.1249999999999998e+00,	6.7949762798203542e-02,	9.7559223176554566e-01],
    [1.2499999999999998e+00,	7.5499736442448381e-02,	9.7559223176554566e-01],
    [1.3750000000000000e+00,	8.3049710086693221e-02,	9.7559223176554566e-01],
    [1.5000000000000000e+00,	9.0599683730938060e-02,	9.7559223176554566e-01],
    [1.6250000000000000e+00,	9.8149657375182886e-02,	9.7559223176554566e-01],
    [1.7499999999999998e+00,	1.0569963101942774e-01,	9.7559223176554566e-01],
    [1.8750000000000000e+00,	1.1324960466367257e-01,	9.7559223176554566e-01],
    [1.9583333333333330e+00,	1.1828292042650249e-01,	9.7559223176554566e-01],
    [2.0000000000000000e+00,	1.2079957830791742e-01,	9.7559223176554566e-01],
    [9.8883384585137646e-01,	1.8487074715047677e-01,	9.3287863735525056e-01],
    [1.0300352560951835e+00,	1.9257369494841331e-01,	9.3287863735525056e-01],
    [1.1124380765827981e+00,	2.0797959054428633e-01,	9.3287863735525056e-01],
    [1.2360423073142206e+00,	2.3108843393809597e-01,	9.3287863735525056e-01],
    [1.3596465380456424e+00,	2.5419727733190556e-01,	9.3287863735525056e-01],
    [1.4832507687770644e+00,	2.7730612072571520e-01,	9.3287863735525056e-01],
    [1.6068549995084866e+00,	3.0041496411952473e-01,	9.3287863735525056e-01],
    [1.7304592302399082e+00,	3.2352380751333432e-01,	9.3287863735525056e-01],
    [1.8540634609713307e+00,	3.4663265090714396e-01,	9.3287863735525056e-01],
    [1.9364662814589453e+00,	3.6203854650301709e-01,	9.3287863735525056e-01],
    [1.9776676917027529e+00,	3.6974149430095354e-01,	9.3287863735525056e-01],
    [9.3541792596869744e-01,	3.7179255260460259e-01,	8.8711407191564851e-01],
    [9.7439367288405998e-01,	3.8728390896312764e-01,	8.8711407191564851e-01],
    [1.0523451667147843e+00,	4.1826662168017775e-01,	8.8711407191564873e-01],
    [1.1692724074608718e+00,	4.6474069075575303e-01,	8.8711407191564873e-01],
    [1.2861996482069591e+00,	5.1121475983132847e-01,	8.8711407191564851e-01],
    [1.4031268889530464e+00,	5.5768882890690386e-01,	8.8711407191564851e-01],
    [1.5200541296991332e+00,	6.0416289798247913e-01,	8.8711407191564851e-01],
    [1.6369813704452205e+00,	6.5063696705805452e-01,	8.8711407191564851e-01],
    [1.7539086111913078e+00,	6.9711103613362968e-01,	8.8711407191564851e-01],
    [1.8318601050220324e+00,	7.2809374885067990e-01,	8.8711407191564851e-01],
    [1.8708358519373949e+00,	7.4358510520920518e-01,	8.8711407191564851e-01],
    [8.4247562770426876e-01,	5.5166140192753399e-01,	8.5965533265188732e-01],
    [8.7757877885861335e-01,	5.7464729367451450e-01,	8.5965533265188732e-01],
    [9.4778508116730209e-01,	6.2061907716847553e-01,	8.5965533265188743e-01],
    [1.0530945346303358e+00,	6.8957675240941751e-01,	8.5965533265188743e-01],
    [1.1584039880933694e+00,	7.5853442765035906e-01,	8.5965533265188732e-01],
    [1.2637134415564031e+00,	8.2749210289130104e-01,	8.5965533265188732e-01],
    [1.3690228950194367e+00,	8.9644977813224269e-01,	8.5965533265188732e-01],
    [1.4743323484824702e+00,	9.6540745337318445e-01,	8.5965533265188732e-01],
    [1.5796418019455041e+00,	1.0343651286141262e+00,	8.5965533265188732e-01],
    [1.6498481042541928e+00,	1.0803369121080875e+00,	8.5965533265188732e-01],
    [1.6849512554085375e+00,	1.1033228038550680e+00,	8.5965533265188732e-01],
    [7.1217992913863115e-01,	7.1217992913863115e-01,	8.5050241956396699e-01],
    [7.4185409285274084e-01,	7.4185409285274084e-01,	8.5050241956396699e-01],
    [8.0120242028095978e-01,	8.0120242028095978e-01,	8.5050241956396710e-01],
    [8.9022491142328886e-01,	8.9022491142328897e-01,	8.5050241956396710e-01],
    [9.7924740256561782e-01,	9.7924740256561771e-01,	8.5050241956396699e-01],
    [1.0682698937079469e+00,	1.0682698937079469e+00,	8.5050241956396699e-01],
    [1.1572923848502756e+00,	1.1572923848502756e+00,	8.5050241956396699e-01],
    [1.2463148759926044e+00,	1.2463148759926046e+00,	8.5050241956396699e-01],
    [1.3353373671349333e+00,	1.3353373671349333e+00,	8.5050241956396699e-01],
    [1.3946856945631525e+00,	1.3946856945631525e+00,	8.5050241956396699e-01],
    [1.4243598582772623e+00,	1.4243598582772623e+00,	8.5050241956396699e-01],
    [5.5166140192753388e-01,	8.4247562770426865e-01,	8.5965533265188743e-01],
    [5.7464729367451450e-01,	8.7757877885861324e-01,	8.5965533265188743e-01],
    [6.2061907716847553e-01,	9.4778508116730209e-01,	8.5965533265188743e-01],
    [6.8957675240941740e-01,	1.0530945346303358e+00,	8.5965533265188743e-01],
    [7.5853442765035894e-01,	1.1584039880933692e+00,	8.5965533265188743e-01],
    [8.2749210289130093e-01,	1.2637134415564031e+00,	8.5965533265188743e-01],
    [8.9644977813224258e-01,	1.3690228950194367e+00,	8.5965533265188743e-01],
    [9.6540745337318434e-01,	1.4743323484824702e+00,	8.5965533265188743e-01],
    [1.0343651286141260e+00,	1.5796418019455036e+00,	8.5965533265188743e-01],
    [1.0803369121080872e+00,	1.6498481042541926e+00,	8.5965533265188743e-01],
    [1.1033228038550678e+00,	1.6849512554085373e+00,	8.5965533265188743e-01],
    [3.7179255260460248e-01,	9.3541792596869744e-01,	8.8711407191564851e-01],
    [3.8728390896312759e-01,	9.7439367288405998e-01,	8.8711407191564851e-01],
    [4.1826662168017781e-01,	1.0523451667147845e+00,	8.8711407191564862e-01],
    [4.6474069075575308e-01,	1.1692724074608716e+00,	8.8711407191564862e-01],
    [5.1121475983132847e-01,	1.2861996482069591e+00,	8.8711407191564851e-01],
    [5.5768882890690386e-01,	1.4031268889530464e+00,	8.8711407191564851e-01],
    [6.0416289798247913e-01,	1.5200541296991334e+00,	8.8711407191564851e-01],
    [6.5063696705805441e-01,	1.6369813704452205e+00,	8.8711407191564851e-01],
    [6.9711103613362968e-01,	1.7539086111913078e+00,	8.8711407191564851e-01],
    [7.2809374885067990e-01,	1.8318601050220324e+00,	8.8711407191564851e-01],
    [7.4358510520920496e-01,	1.8708358519373949e+00,	8.8711407191564851e-01],
    [1.8487074715047677e-01,	9.8883384585137646e-01,	9.3287863735525045e-01],
    [1.9257369494841331e-01,	1.0300352560951838e+00,	9.3287863735525045e-01],
    [2.0797959054428636e-01,	1.1124380765827984e+00,	9.3287863735525045e-01],
    [2.3108843393809597e-01,	1.2360423073142206e+00,	9.3287863735525045e-01],
    [2.5419727733190556e-01,	1.3596465380456426e+00,	9.3287863735525045e-01],
    [2.7730612072571514e-01,	1.4832507687770649e+00,	9.3287863735525045e-01],
    [3.0041496411952479e-01,	1.6068549995084866e+00,	9.3287863735525045e-01],
    [3.2352380751333432e-01,	1.7304592302399087e+00,	9.3287863735525045e-01],
    [3.4663265090714401e-01,	1.8540634609713309e+00,	9.3287863735525045e-01],
    [3.6203854650301698e-01,	1.9364662814589455e+00,	9.3287863735525045e-01],
    [3.6974149430095354e-01,	1.9776676917027529e+00,	9.3287863735525045e-01],
    [6.0399789153958702e-02,	1.0000000000000000e+00,	9.7559223176554566e-01],
    [6.2916447035373635e-02,	1.0416666666666667e+00,	9.7559223176554566e-01],
    [6.7949762798203542e-02,	1.1249999999999998e+00,	9.7559223176554566e-01],
    [7.5499736442448367e-02,	1.2499999999999998e+00,	9.7559223176554566e-01],
    [8.3049710086693221e-02,	1.3750000000000000e+00,	9.7559223176554566e-01],
    [9.0599683730938046e-02,	1.5000000000000000e+00,	9.7559223176554566e-01],
    [9.8149657375182886e-02,	1.6250000000000000e+00,	9.7559223176554566e-01],
    [1.0569963101942773e-01,	1.7499999999999998e+00,	9.7559223176554566e-01],
    [1.1324960466367257e-01,	1.8750000000000000e+00,	9.7559223176554566e-01],
    [1.1828292042650244e-01,	1.9583333333333330e+00,	9.7559223176554566e-01],
    [1.2079957830791740e-01,	2.0000000000000000e+00,	9.7559223176554566e-01],
    [0.0000000000000000e+00,	1.0000000000000000e+00,	1.0000000000000000e+00],
    [0.0000000000000000e+00,	1.0416666666666667e+00,	1.0000000000000000e+00],
    [0.0000000000000000e+00,	1.1249999999999998e+00,	1.0000000000000000e+00],
    [0.0000000000000000e+00,	1.2500000000000000e+00,	1.0000000000000000e+00],
    [0.0000000000000000e+00,	1.3750000000000000e+00,	1.0000000000000000e+00],
    [0.0000000000000000e+00,	1.5000000000000000e+00,	1.0000000000000000e+00],
    [0.0000000000000000e+00,	1.6250000000000000e+00,	1.0000000000000000e+00],
    [0.0000000000000000e+00,	1.7500000000000000e+00,	1.0000000000000000e+00],
    [0.0000000000000000e+00,	1.8750000000000000e+00,	1.0000000000000000e+00],
    [0.0000000000000000e+00,	1.9583333333333333e+00,	1.0000000000000000e+00],
    [0.0000000000000000e+00,	2.0000000000000000e+00,	1.0000000000000000e+00]])

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

# 座標，パッチコネクティビティー，境界の辺をtxt出力
# [patch_number, xi_or_eta, 0_or_1(int)] (auto marge)
boundary_array_0 = np.array([[0, eta, 0],
                             [0, xi, 1]])

boundary_array_1 = np.array([[0,  xi, 1],
                             [0, eta, 1]])

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
ax1.set_xlim(-1, 3)
ax1.set_ylim(-1, 3)
plt.show()