import math
import numpy as np

division_ele_xi = 10.
calc_xi_loc = np.zeros((1000))

knot_vec_xi_loc = np.array([0.000000000000000000000e+00,   0.000000000000000000000e+00,   0.000000000000000000000e+00,   6.250000000000000000000e-02,   1.250000000000000000000e-01,   1.875000000000000000000e-01,   2.500000000000000000000e-01,   3.125000000000000000000e-01,   3.750000000000000000000e-01,   4.375000000000000000000e-01,   5.000000000000000000000e-01,   5.625000000000000000000e-01,   6.250000000000000000000e-01,   6.875000000000000000000e-01,   7.500000000000000000000e-01,   8.125000000000000000000e-01,   8.750000000000000000000e-01,   9.375000000000000000000e-01,   1.000000000000000000000e+00,   1.000000000000000000000e+00,   1.000000000000000000000e+00])
knot_n_xi_loc = knot_vec_xi_loc.shape[0]

# print(knot_vec_xi_loc)
# print(knot_n_xi_loc)

k = 0
l = 0

for i in range(knot_n_xi_loc-1):
    if(knot_vec_xi_loc[i] != knot_vec_xi_loc[i+1]):
        calc_xi_loc[k] = knot_vec_xi_loc[i]
        print(calc_xi_loc[k])
        k += 1
        l += 1
        if(int(division_ele_xi) > 1):
            temp1 = (knot_vec_xi_loc[i + 1] - knot_vec_xi_loc[i]) / float(division_ele_xi)
            for j in range(int(division_ele_xi-1)):
                calc_xi_loc[k] = calc_xi_loc[k-1] + temp1
                k += 1
calc_xi_loc[k] = knot_vec_xi_loc[knot_n_xi_loc - 1]
division_n_xi = k + 1
element_n_xi = l
print(division_n_xi)
print(element_n_xi)


#-----------------------------------------------------------

division_ele_eta = 10.
calc_eta_loc = np.zeros((1000))

knot_vec_eta_loc = np.array([0.000000000000000000000e+00,   0.000000000000000000000e+00,   0.000000000000000000000e+00,   1.250000000000000000000e-01,   2.500000000000000000000e-01,   3.750000000000000000000e-01,   5.000000000000000000000e-01,   6.250000000000000000000e-01,   7.500000000000000000000e-01,   8.750000000000000000000e-01,   1.000000000000000000000e+00,   1.000000000000000000000e+00,   1.000000000000000000000e+00])
knot_n_eta_loc = knot_vec_eta_loc.shape[0]

# print(knot_vec_eta_loc)
# print(knot_n_eta_loc)

k = 0
l = 0

for i in range(knot_n_eta_loc-1):
    if(knot_vec_eta_loc[i] != knot_vec_eta_loc[i+1]):
        calc_eta_loc[k] = knot_vec_eta_loc[i]
        print(calc_eta_loc[k])
        k += 1
        l += 1
        if(int(division_ele_eta) > 1):
            temp1 = (knot_vec_eta_loc[i + 1] - knot_vec_eta_loc[i]) / float(division_ele_eta)
            for j in range(int(division_ele_eta-1)):
                calc_eta_loc[k] = calc_eta_loc[k-1] + temp1
                k += 1
calc_eta_loc[k] = knot_vec_eta_loc[knot_n_eta_loc - 1]
division_n_eta = k + 1
element_n_eta = l
print(division_n_eta)
print(element_n_eta)


#-----------------------------------------------------------

for i in range(int(division_n_xi)):
    ii = int(i // division_ele_xi)
    kk = int(i % division_ele_xi)
    for j in range(int(division_n_eta)):
        jj = int(j // division_ele_eta)
        ll = int(j % division_ele_eta)