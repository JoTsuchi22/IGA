import numpy as np
import matplotlib.pyplot as plt
import math
import make_inputdate_function as mif

# difine filename
filename = "auto_test"

# patch info (ξ方向のコントロールポイント個数, η方向のコントロールポイント個数)
patch_info = np.array([[6, 6],
                       [6, 6],
                       [6, 6],
                       [6, 6],
                       [6, 6]])

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
   [[1.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000],
    [1.25000000000000000000, 0.00000000000000000000, 1.00000000000000000000],
    [1.75000000000000000000, 0.00000000000000000000, 1.00000000000000000000],
    [2.25000000000000000000, 0.00000000000000000000, 1.00000000000000000000],
    [2.75000000000000000000, 0.00000000000000000000, 1.00000000000000000000],
    [3.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000],
    [1.00000000000000000000, 0.09752680458061159519, 0.98096988312782174013],
    [1.25423332871344195283, 0.13279827163072491225, 0.98334864773684405037],
    [1.75902774164546205071, 0.20283173102377569097, 0.98810617695488855983],
    [2.25898448319580280597, 0.27219402834849487638, 0.99286370617293318031],
    [2.75417276408932698573, 0.34089476566537479929, 0.99762123539097768976],
    [3.00000000000000000000, 0.37500000000000000000, 1.00000000000000000000],
    [0.96193976625564348026, 0.29079789987237392168, 0.96193976625564348026],
    [1.22547368779538601302, 0.39866544331559600201, 0.96669729547368787870],
    [1.74483560299975493102, 0.61124640300284283079, 0.97621235390977711965],
    [2.25417089439158591091, 0.81972334726251783188, 0.98572741234586636061],
    [2.75376714187515680976, 1.02421398594250900871, 0.99524247078195537952],
    [3.00000000000000000000, 1.12500000000000000000, 1.00000000000000000000],
    [0.88581929876693032977, 0.47456896485780580841, 0.96193976625564336924],
    [1.15919607233557697690, 0.65565344726240026851, 0.96669729547368798972],
    [1.69795587928358493457, 1.01252736712041291334, 0.97621235390977711965],
    [2.22631457334395088665, 1.36251160092681389280, 0.98572741234586636061],
    [2.74457047537766607448, 1.70580375610083700977, 0.99524247078195537952],
    [3.00000000000000000000, 1.87500000000000000000, 1.00000000000000000000],
    [0.77606864605295333703, 0.63814491632014191946, 0.98096988312782174013],
    [1.05876738152502802848, 0.89070730711771917143, 0.98334864773684405037],
    [1.62008144263343289992, 1.39218397989130493464, 0.98810617695488855983],
    [2.17601617920893497526, 1.88885477667596091145, 0.99286370617293318031],
    [2.72664855120205817940, 2.38078845331195498858, 0.99762123539097768976],
    [3.00000000000000000000, 2.62500000000000000000, 1.00000000000000000000],
    [0.70710678118654757274, 0.70710678118654757274, 1.00000000000000000000],
    [0.99371843353822875144, 0.99371843353822897349, 1.00000000000000000000],
    [1.56694173824159199704, 1.56694173824159199704, 1.00000000000000000000],
    [2.14016504294495479854, 2.14016504294495479854, 1.00000000000000000000],
    [2.71338834764831782209, 2.71338834764831782209, 1.00000000000000000000],
    [3.00000000000000000000, 3.00000000000000000000, 1.00000000000000000000]])

num = 1
patch[num,:patch_info[num][0]*patch_info[num][1],:] = np.array(
   [[3.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000],
    [3.87500000000000000000, 0.00000000000000000000, 1.00000000000000000000],
    [5.62500000000000000000, 0.00000000000000000000, 1.00000000000000000000],
    [7.37500000000000000000, 0.00000000000000000000, 1.00000000000000000000],
    [9.12500000000000000000, 0.00000000000000000000, 1.00000000000000000000],
    [10.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000],
    [3.00000000000000000000, 0.37500000000000000000, 1.00000000000000000000],
    [3.87500000000000000000, 0.37500000000000000000, 1.00000000000000000000],
    [5.62500000000000000000, 0.37500000000000000000, 1.00000000000000000000],
    [7.37500000000000000000, 0.37500000000000000000, 1.00000000000000000000],
    [9.12500000000000000000, 0.37500000000000000000, 1.00000000000000000000],
    [10.00000000000000000000, 0.37500000000000000000, 1.00000000000000000000],
    [3.00000000000000000000, 1.12500000000000000000, 1.00000000000000000000],
    [3.87500000000000000000, 1.12500000000000000000, 1.00000000000000000000],
    [5.62500000000000000000, 1.12500000000000000000, 1.00000000000000000000],
    [7.37500000000000000000, 1.12500000000000000000, 1.00000000000000000000],
    [9.12500000000000000000, 1.12500000000000000000, 1.00000000000000000000],
    [10.00000000000000000000, 1.12500000000000000000, 1.00000000000000000000],
    [3.00000000000000000000, 1.87500000000000000000, 1.00000000000000000000],
    [3.87500000000000000000, 1.87500000000000000000, 1.00000000000000000000],
    [5.62500000000000000000, 1.87500000000000000000, 1.00000000000000000000],
    [7.37500000000000000000, 1.87500000000000000000, 1.00000000000000000000],
    [9.12500000000000000000, 1.87500000000000000000, 1.00000000000000000000],
    [10.00000000000000000000, 1.87500000000000000000, 1.00000000000000000000],
    [3.00000000000000000000, 2.62500000000000000000, 1.00000000000000000000],
    [3.87500000000000000000, 2.62500000000000000000, 1.00000000000000000000],
    [5.62500000000000000000, 2.62500000000000000000, 1.00000000000000000000],
    [7.37500000000000000000, 2.62500000000000000000, 1.00000000000000000000],
    [9.12500000000000000000, 2.62500000000000000000, 1.00000000000000000000],
    [10.00000000000000000000, 2.62500000000000000000, 1.00000000000000000000],
    [3.00000000000000000000, 3.00000000000000000000, 1.00000000000000000000],
    [3.87500000000000000000, 3.00000000000000000000, 1.00000000000000000000],
    [5.62500000000000000000, 3.00000000000000000000, 1.00000000000000000000],
    [7.37500000000000000000, 3.00000000000000000000, 1.00000000000000000000],
    [9.12500000000000000000, 3.00000000000000000000, 1.00000000000000000000],
    [10.00000000000000000000, 3.00000000000000000000, 1.00000000000000000000]])

num = 2
patch[num,:patch_info[num][0]*patch_info[num][1],:] = np.array(
   [[3.00000000000000000000, 3.00000000000000000000, 1.00000000000000000000],
    [3.87500000000000000000, 3.00000000000000000000, 1.00000000000000000000],
    [5.62500000000000000000, 3.00000000000000000000, 1.00000000000000000000],
    [7.37500000000000000000, 3.00000000000000000000, 1.00000000000000000000],
    [9.12500000000000000000, 3.00000000000000000000, 1.00000000000000000000],
    [10.00000000000000000000, 3.00000000000000000000, 1.00000000000000000000],
    [3.00000000000000000000, 3.87500000000000000000, 1.00000000000000000000],
    [3.87500000000000000000, 3.87500000000000000000, 1.00000000000000000000],
    [5.62500000000000000000, 3.87500000000000000000, 1.00000000000000000000],
    [7.37500000000000000000, 3.87500000000000000000, 1.00000000000000000000],
    [9.12500000000000000000, 3.87500000000000000000, 1.00000000000000000000],
    [10.00000000000000000000, 3.87500000000000000000, 1.00000000000000000000],
    [3.00000000000000000000, 5.62500000000000000000, 1.00000000000000000000],
    [3.87500000000000000000, 5.62500000000000000000, 1.00000000000000000000],
    [5.62500000000000000000, 5.62500000000000000000, 1.00000000000000000000],
    [7.37500000000000000000, 5.62500000000000000000, 1.00000000000000000000],
    [9.12500000000000000000, 5.62500000000000000000, 1.00000000000000000000],
    [10.00000000000000000000, 5.62500000000000000000, 1.00000000000000000000],
    [3.00000000000000000000, 7.37500000000000000000, 1.00000000000000000000],
    [3.87500000000000000000, 7.37500000000000000000, 1.00000000000000000000],
    [5.62500000000000000000, 7.37500000000000000000, 1.00000000000000000000],
    [7.37500000000000000000, 7.37500000000000000000, 1.00000000000000000000],
    [9.12500000000000000000, 7.37500000000000000000, 1.00000000000000000000],
    [10.00000000000000000000, 7.37500000000000000000, 1.00000000000000000000],
    [3.00000000000000000000, 9.12500000000000000000, 1.00000000000000000000],
    [3.87500000000000000000, 9.12500000000000000000, 1.00000000000000000000],
    [5.62500000000000000000, 9.12500000000000000000, 1.00000000000000000000],
    [7.37500000000000000000, 9.12500000000000000000, 1.00000000000000000000],
    [9.12500000000000000000, 9.12500000000000000000, 1.00000000000000000000],
    [10.00000000000000000000, 9.12500000000000000000, 1.00000000000000000000],
    [3.00000000000000000000, 10.00000000000000000000, 1.00000000000000000000],
    [3.87500000000000000000, 10.00000000000000000000, 1.00000000000000000000],
    [5.62500000000000000000, 10.00000000000000000000, 1.00000000000000000000],
    [7.37500000000000000000, 10.00000000000000000000, 1.00000000000000000000],
    [9.12500000000000000000, 10.00000000000000000000, 1.00000000000000000000],
    [10.00000000000000000000, 10.00000000000000000000, 1.00000000000000000000]])

num = 3
patch[num,:patch_info[num][0]*patch_info[num][1],:] = np.array(
   [[0.00000000000000000000, 3.00000000000000000000, 1.00000000000000000000],
    [0.37500000000000000000, 3.00000000000000000000, 1.00000000000000000000],
    [1.12500000000000000000, 3.00000000000000000000, 1.00000000000000000000],
    [1.87500000000000000000, 3.00000000000000000000, 1.00000000000000000000],
    [2.62500000000000000000, 3.00000000000000000000, 1.00000000000000000000],
    [3.00000000000000000000, 3.00000000000000000000, 1.00000000000000000000],
    [0.00000000000000000000, 3.87500000000000000000, 1.00000000000000000000],
    [0.37500000000000000000, 3.87500000000000000000, 1.00000000000000000000],
    [1.12500000000000000000, 3.87500000000000000000, 1.00000000000000000000],
    [1.87500000000000000000, 3.87500000000000000000, 1.00000000000000000000],
    [2.62500000000000000000, 3.87500000000000000000, 1.00000000000000000000],
    [3.00000000000000000000, 3.87500000000000000000, 1.00000000000000000000],
    [0.00000000000000000000, 5.62500000000000000000, 1.00000000000000000000],
    [0.37500000000000000000, 5.62500000000000000000, 1.00000000000000000000],
    [1.12500000000000000000, 5.62500000000000000000, 1.00000000000000000000],
    [1.87500000000000000000, 5.62500000000000000000, 1.00000000000000000000],
    [2.62500000000000000000, 5.62500000000000000000, 1.00000000000000000000],
    [3.00000000000000000000, 5.62500000000000000000, 1.00000000000000000000],
    [0.00000000000000000000, 7.37500000000000000000, 1.00000000000000000000],
    [0.37500000000000000000, 7.37500000000000000000, 1.00000000000000000000],
    [1.12500000000000000000, 7.37500000000000000000, 1.00000000000000000000],
    [1.87500000000000000000, 7.37500000000000000000, 1.00000000000000000000],
    [2.62500000000000000000, 7.37500000000000000000, 1.00000000000000000000],
    [3.00000000000000000000, 7.37500000000000000000, 1.00000000000000000000],
    [0.00000000000000000000, 9.12500000000000000000, 1.00000000000000000000],
    [0.37500000000000000000, 9.12500000000000000000, 1.00000000000000000000],
    [1.12500000000000000000, 9.12500000000000000000, 1.00000000000000000000],
    [1.87500000000000000000, 9.12500000000000000000, 1.00000000000000000000],
    [2.62500000000000000000, 9.12500000000000000000, 1.00000000000000000000],
    [3.00000000000000000000, 9.12500000000000000000, 1.00000000000000000000],
    [0.00000000000000000000, 10.00000000000000000000, 1.00000000000000000000],
    [0.37500000000000000000, 10.00000000000000000000, 1.00000000000000000000],
    [1.12500000000000000000, 10.00000000000000000000, 1.00000000000000000000],
    [1.87500000000000000000, 10.00000000000000000000, 1.00000000000000000000],
    [2.62500000000000000000, 10.00000000000000000000, 1.00000000000000000000],
    [3.00000000000000000000, 10.00000000000000000000, 1.00000000000000000000]])

num = 4
patch[num,:patch_info[num][0]*patch_info[num][1],:] = np.array(
   [[0.70710678118654768376, 0.70710678118654757274, 1.00000000000000000000],
    [0.99371843353822908451, 0.99371843353822897349, 1.00000000000000000000],
    [1.56694173824159199704, 1.56694173824159199704, 1.00000000000000000000],
    [2.14016504294495479854, 2.14016504294495479854, 1.00000000000000000000],
    [2.71338834764831782209, 2.71338834764831782209, 1.00000000000000000000],
    [3.00000000000000000000, 3.00000000000000000000, 1.00000000000000000000],
    [0.63814491632014214151, 0.77606864605295333703, 0.98096988312782174013],
    [0.89070730711771939347, 1.05876738152502802848, 0.98334864773684405037],
    [1.39218397989130604486, 1.62008144263343289992, 0.98810617695488855983],
    [1.88885477667596091145, 2.17601617920893497526, 0.99286370617293318031],
    [2.38078845331195498858, 2.72664855120205817940, 0.99762123539097768976],
    [2.62500000000000000000, 3.00000000000000000000, 1.00000000000000000000],
    [0.47456896485780608597, 0.88581929876693032977, 0.96193976625564336924],
    [0.65565344726240060158, 1.15919607233557697690, 0.96669729547368798972],
    [1.01252736712041291334, 1.69795587928358493457, 0.97621235390977711965],
    [1.36251160092681389280, 2.22631457334395088665, 0.98572741234586636061],
    [1.70580375610083700977, 2.74457047537766607448, 0.99524247078195537952],
    [1.87500000000000000000, 3.00000000000000000000, 1.00000000000000000000],
    [0.29079789987237431026, 0.96193976625564348026, 0.96193976625564348026],
    [0.39866544331559622405, 1.22547368779538601302, 0.96669729547368787870],
    [0.61124640300284294181, 1.74483560299975493102, 0.97621235390977711965],
    [0.81972334726251783188, 2.25417089439158591091, 0.98572741234586636061],
    [1.02421398594250900871, 2.75376714187515680976, 0.99524247078195537952],
    [1.12500000000000000000, 3.00000000000000000000, 1.00000000000000000000],
    [0.09752680458061180335, 1.00000000000000000000, 0.98096988312782174013],
    [0.13279827163072510654, 1.25423332871344195283, 0.98334864773684405037],
    [0.20283173102377580199, 1.75902774164546205071, 0.98810617695488855983],
    [0.27219402834849487638, 2.25898448319580280597, 0.99286370617293318031],
    [0.34089476566537479929, 2.75417276408932698573, 0.99762123539097768976],
    [0.37500000000000000000, 3.00000000000000000000, 1.00000000000000000000],
    [0.00000000000000012246, 1.00000000000000000000, 1.00000000000000000000],
    [0.00000000000000010716, 1.25000000000000000000, 1.00000000000000000000],
    [0.00000000000000007654, 1.75000000000000000000, 1.00000000000000000000],
    [0.00000000000000004592, 2.25000000000000000000, 1.00000000000000000000],
    [0.00000000000000001531, 2.75000000000000000000, 1.00000000000000000000],
    [0.00000000000000000000, 3.00000000000000000000, 1.00000000000000000000]])


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
connection_vec = np.array([patch_number, xi, 0, connect_patch_number, xi, 1, positive])
globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool = \
    mif.connect_patch_arg_1boundary(patch_number, connection_vec, patch_info, patch, patch_bool, globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool)

patch_number = 2
connect_patch_number = 1
connection_vec = np.array([patch_number, eta, 0, connect_patch_number, eta, 1, positive])
globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool = \
    mif.connect_patch_arg_1boundary(patch_number, connection_vec, patch_info, patch, patch_bool, globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool)

patch_number = 3
connect_patch_number = 2
connection_vec = np.array([patch_number, xi, 1, connect_patch_number, xi, 0, positive])
globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool = \
    mif.connect_patch_arg_1boundary(patch_number, connection_vec, patch_info, patch, patch_bool, globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool)

patch_number = 4
connect_patch_number_1 = 3
connect_patch_number_2 = 0
connection_vec = np.array([[patch_number, xi, 1, connect_patch_number_1, eta, 0, negative],
                           [patch_number, eta, 0, connect_patch_number_2, eta, 1, positive]])
globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool = \
    mif.connect_patch_arg_2boundary(patch_number, connection_vec, patch_info, patch, patch_bool, globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool)


# 座標，パッチコネクティビティー，境界の辺をtxt出力
# [patch_number, xi_or_eta, 0_or_1(int)] (auto marge)
boundary_array_0 = np.array([[0, eta, 0],
                             [1, eta, 0]])

boundary_array_1 = np.array([[3,  xi, 0],
                             [4, eta, 1]])

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