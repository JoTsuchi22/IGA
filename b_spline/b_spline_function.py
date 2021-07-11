import numpy as np
import matplotlib.pyplot as plt
import math


def def_knot(m, n):
    knot = np.zeros((m))
    for i in range(m-2*(n+1)):
        knot[i+n+1] = (i + 1) / (m - 2 * (n + 1) + 1)
    for i in range(n+1):
        knot[m-i-1] = 1.0
    return knot


def def_knot_C0(m, n, make_C0_CP):
    knot = np.zeros((2*(n+1)))
    knot_var = np.zeros((m-2*(n+1)-make_C0_CP.shape[0]))
    for i in range(n+1):
        knot[2*(n+1)-i-1] = 1.0
    for i in range(m-2*(n+1)-make_C0_CP.shape[0]):
        knot_var[i] = (i + 1) / (m - 2 * (n + 1) + 1 - make_C0_CP.shape[0])
    for i in range(make_C0_CP.shape[0]):
        knot_var = np.append(
            knot_var, knot_var[make_C0_CP[i]-make_C0_CP.shape[0]-1])
    knot = np.sort(np.append(knot, knot_var))
    return knot


def knot_insertion_A(CP, n, l, knot, new_knot_position):
    knot_rate = 0.5
    number_of_insertion = new_knot_position.shape[0]
    new_knot = np.zeros((number_of_insertion))
    knot_var_old = knot
    for i in range(number_of_insertion):
        new_knot[i] = knot_rate * (knot[n+new_knot_position[i]] -
                                   knot[n+new_knot_position[i]-1]) + knot[n+new_knot_position[i]-1]
        new_knot_number = n + new_knot_position[i] + i
        knot_var_old = np.insert(knot_var_old, new_knot_number, new_knot[i])
    knot_var = knot_var_old

    T = np.zeros((n+1, l+number_of_insertion, l))
    for i in range(l+number_of_insertion):
        for k in range(n+1):
            for j in range(l):
                if k == 0:
                    if (knot[j] <= knot_var[i]) and (knot_var[i] < knot[j+1]):
                        T[k][i][j] = 1.0
                    else:
                        T[k][i][j] = 0.0
                else:
                    if (knot[j+k] - knot[j]) == 0:
                        a = 0.0
                    if (knot[j+k] - knot[j]) != 0:
                        a = ((knot_var[i+k] - knot[j]) /
                             (knot[j+k] - knot[j])) * T[k-1][i][j]
                    if (knot[j+k+1] - knot[j+1]) == 0:
                        b = 0.0
                    if (knot[j+k+1] - knot[j+1]) != 0:
                        b = ((knot[j+k+1] - knot_var[i+k]) /
                             (knot[j+k+1] - knot[j+1])) * T[k-1][i][j+1]
                    T[k][i][j] = a + b
    CP = np.array(np.dot(T[n, :, :], CP))
    knot = np.array(knot_var)
    l = l + number_of_insertion
    m = l + n + 1
    return CP, l, m, knot


def knot_insertion_B(CP, n, l, knot, insert_knot):
    number_of_insertion = insert_knot.shape[0]
    knot_var_old = knot
    knot_var_old = np.append(knot_var_old, insert_knot)
    knot_var = np.sort(knot_var_old)

    T = np.zeros((n+1, l+number_of_insertion, l))
    for i in range(l+number_of_insertion):
        for k in range(n+1):
            for j in range(l):
                if k == 0:
                    if (knot[j] <= knot_var[i]) and (knot_var[i] < knot[j+1]):
                        T[k][i][j] = 1.0
                    else:
                        T[k][i][j] = 0.0
                else:
                    if (knot[j+k] - knot[j]) == 0:
                        a = 0.0
                    if (knot[j+k] - knot[j]) != 0:
                        a = ((knot_var[i+k] - knot[j]) /
                             (knot[j+k] - knot[j])) * T[k-1][i][j]
                    if (knot[j+k+1] - knot[j+1]) == 0:
                        b = 0.0
                    if (knot[j+k+1] - knot[j+1]) != 0:
                        b = ((knot[j+k+1] - knot_var[i+k]) /
                             (knot[j+k+1] - knot[j+1])) * T[k-1][i][j+1]
                    T[k][i][j] = a + b
    CP = np.array(np.dot(T[n, :, :], CP))
    knot = np.array(knot_var)
    l = l + number_of_insertion
    m = l + n + 1
    return CP, l, m, knot


def knot_removal(CP, n, l, knot, removal_knot):
    knot_bool_vec_1 = np.ones((knot.shape[0]), dtype=bool)
    knot_bool_vec_2 = np.zeros(
        (knot.shape[0], np.unique(removal_knot).shape[0]), dtype=bool)
    length_vec = np.zeros(np.unique(removal_knot).shape[0])
    for i in range(np.unique(removal_knot).shape[0]):
        knot_bool_vec_3 = np.zeros((removal_knot.shape[0]), dtype=bool)
        for j in range(removal_knot.shape[0]):
            if np.unique(removal_knot[i]) == removal_knot[j]:
                knot_bool_vec_3[j] = True
        length_vec[i] = removal_knot[knot_bool_vec_3].shape[0]
        for j in range(knot.shape[0]):
            if removal_knot[i] == knot[j]:
                knot_bool_vec_1[j] = False
                knot_bool_vec_2[j][i] = True
    knot_var = knot[knot_bool_vec_1]
    for i in range(np.unique(removal_knot).shape[0]):
        for j in range(int(knot[knot_bool_vec_2[:, i]].shape[0]-length_vec[i])):
            knot_var = np.append(knot_var, np.unique(removal_knot)[i])
    knot_var_old = np.sort(knot_var)
    number_of_removal = removal_knot.shape[0]
    knot_var = knot
    knot = knot_var_old

    l = l - number_of_removal
    T = np.zeros((n+1, l+number_of_removal, l))
    for i in range(l+number_of_removal):
        for k in range(n+1):
            for j in range(l):
                if k == 0:
                    if (knot[j] <= knot_var[i]) and (knot_var[i] < knot[j+1]):
                        T[k][i][j] = 1.0
                    else:
                        T[k][i][j] = 0.0
                else:
                    if (knot[j+k] - knot[j]) == 0:
                        a = 0.0
                    if (knot[j+k] - knot[j]) != 0:
                        a = ((knot_var[i+k] - knot[j]) /
                             (knot[j+k] - knot[j])) * T[k-1][i][j]
                    if (knot[j+k+1] - knot[j+1]) == 0:
                        b = 0.0
                    if (knot[j+k+1] - knot[j+1]) != 0:
                        b = ((knot[j+k+1] - knot_var[i+k]) /
                             (knot[j+k+1] - knot[j+1])) * T[k-1][i][j+1]
                    T[k][i][j] = a + b
    CP = np.array(np.dot(np.linalg.pinv(T[n, :, :]), CP))
    m = l + n + 1
    return CP, l, m, knot


def Bézier_order_elevarion_function(CP_Bézier_matrix, n, elevation_degree):
    for j in range(elevation_degree):
        alpha = np.zeros((n))
        new_CP = np.zeros((n+1, 2))
        for i in range(n+1):
            if i != n:
                alpha[i] = (i + 1) / (n + 1)
                a_x = (1 - alpha[i]) * CP_Bézier_matrix[i+1][0]
                a_y = (1 - alpha[i]) * CP_Bézier_matrix[i+1][1]
                b_x = alpha[i] * CP_Bézier_matrix[i][0]
                b_y = alpha[i] * CP_Bézier_matrix[i][1]
                new_CP[i][0] = a_x + b_x
                new_CP[i][1] = a_y + b_y
            if i == n:
                new_CP[i][0] = CP_Bézier_matrix[i][0]
                new_CP[i][1] = CP_Bézier_matrix[i][1]
        if j != elevation_degree-1 and elevation_degree != 1:
            next_CP = np.zeros((n+2, 2))
            next_CP[0][0] = CP_Bézier_matrix[0][0]
            next_CP[0][1] = CP_Bézier_matrix[0][1]
            for i in range(n+1):
                next_CP[i+1][0] = new_CP[i][0]
                next_CP[i+1][1] = new_CP[i][1]
            CP_Bézier_matrix = next_CP
            n = n + 1
    return new_CP


def order_elevation(CP, n, l, knot, elevation_degree):
    if elevation_degree == 0:
        m = l + n + 1
        return CP, l, m, n, knot
    if elevation_degree != 0:
        knot_bool_vec = np.zeros((knot.shape[0]), dtype=bool)
        knot_var_1 = np.zeros((1))
        for i in range(np.unique(knot).shape[0]):
            for j in range(knot.shape[0]):
                if np.unique(knot)[i] == knot[j]:
                    knot_bool_vec[j] = True
                else:
                    knot_bool_vec[j] = False
            if (knot[knot_bool_vec].shape[0] < n):
                knot_var_1 = np.append(knot_var_1, np.unique(knot)[i])
        knot_var_2 = knot_var_1[1:]
        knot_var_3 = knot_var_2
        for i in range(n-2):
            knot_var_2 = np.append(knot_var_2, knot_var_2)
        insert_knot = np.sort(knot_var_2)
        CP, l, m, knot = knot_insertion_B(CP, n, l, knot, insert_knot)
        a = 1 + ((l - 1 - n) // n)
        CP_Bézier_number = np.zeros((a))
        for i in range(a):
            CP_Bézier_number[i] = i * n
        CP_Bézier_matrix = np.zeros((CP_Bézier_number.shape[0], n+1, 2))
        new_CP = CP[0, :]
        for s in range(a):
            for t in range(n+1):
                point_number = int(CP_Bézier_number[s] + t)
                CP_Bézier_matrix[s][t][0] = CP[point_number][0]
                CP_Bézier_matrix[s][t][1] = CP[point_number][1]
        for s in range(a):
            new_CP = np.append(new_CP, Bézier_order_elevarion_function(
                CP_Bézier_matrix[s, :, :], n, elevation_degree))
        CP = np.reshape(new_CP, [l+elevation_degree*a, 2])
        removal_knot = knot_var_3
        for i in range(elevation_degree):
            knot = np.sort(np.append(knot, np.unique(knot)))
        l = l + elevation_degree*a
        n = n + elevation_degree
        CP, l, m, knot = knot_removal(CP, n, l, knot, removal_knot)
        return CP, l, m, n, knot


def basisfunction_return_N(N, delta, knot, l, m, n):
    for i in range(delta):
        nowknot = (i / (delta - 1)) * knot[m-1]
        for k in range(n+1):
            for j in range(l):
                if k == 0:
                    if (knot[j] <= nowknot) and (nowknot < knot[j+1]):
                        N[k][i][j] = 1.0
                    else:
                        N[k][i][j] = 0.0
                    if i == delta-1:
                        if (knot[j] <= nowknot) and (nowknot <= knot[j+1]):
                            N[k][i][j] = 1.0
                        else:
                            N[k][i][j] = 0.0
                else:
                    if (knot[j+k] - knot[j]) == 0:
                        a = 0
                    if (knot[j+k] - knot[j]) != 0:
                        a = ((nowknot - knot[j]) /
                             (knot[j+k] - knot[j])) * N[k-1][i][j]
                    if (knot[j+k+1] - knot[j+1]) == 0:
                        b = 0
                    if (knot[j+k+1] - knot[j+1]) != 0:
                        b = ((knot[j+k+1] - nowknot) /
                             (knot[j+k+1] - knot[j+1])) * N[k-1][i][j+1]
                    N[k][i][j] = a + b
    return N


def basisfunction_return_N_at_xi(N_xi, xi, knot, l, n):
    for i in range(xi.shape[0]):
        nowknot = xi[i]
        for k in range(n+1):
            for j in range(l):
                if k == 0:
                    if (knot[j] <= nowknot) and (nowknot < knot[j+1]):
                        N_xi[k][i][j] = 1.0
                    else:
                        N_xi[k][i][j] = 0.0
                    if i == xi.shape[0]-1:
                        if (knot[j] <= nowknot) and (nowknot <= knot[j+1]):
                            N_xi[k][i][j] = 1.0
                        else:
                            N_xi[k][i][j] = 0.0
                else:
                    if (knot[j+k] - knot[j]) == 0:
                        a = 0
                    if (knot[j+k] - knot[j]) != 0:
                        a = ((nowknot - knot[j]) /
                             (knot[j+k] - knot[j])) * N_xi[k-1][i][j]
                    if (knot[j+k+1] - knot[j+1]) == 0:
                        b = 0
                    if (knot[j+k+1] - knot[j+1]) != 0:
                        b = ((knot[j+k+1] - nowknot) /
                             (knot[j+k+1] - knot[j+1])) * N_xi[k-1][i][j+1]
                    N_xi[k][i][j] = a + b
    return N_xi


def weight_basisfunction_1parameter_return_R(R, N, w, delta, n, l):
    a = np.zeros((delta))
    for i in range(delta):
        for p in range(l):
            a[i] += N[n][i][p] * w[p]
        for p in range(l):
            R[i][p] = (N[n][i][p] * w[p]) / a[i]
    return R


def weight_basisfunction_2parameter_return_R(R, N, M, w, delta, n, l):
    a = np.zeros((delta[0], delta[1]))
    for i in range(delta[0]):
        for j in range(delta[1]):
            for p in range(l[0]):
                for q in range(l[1]):
                    a[i][j] += N[n[0]][i][p] * M[n[1]][j][q] * w[p][q]
            for p in range(l[0]):
                for q in range(l[1]):
                    R[i][j][p][q] = (N[n[0]][i][p] * M[n[1]][j][q] * w[p][q]) / a[i][j]
    return R


def weight_basisfunction_3parameter_return_R(R, N, M, L, w, delta, n, l):
    a = np.zeros((delta[0], delta[1], delta[2]))
    for i in range(delta[0]):
        for j in range(delta[1]):
            for k in range(delta[2]):
                for p in range(l[0]):
                    for q in range(l[1]):
                        for r in range(l[2]):
                            a[i][j][k] += N[n[0]][i][p] * \
                                M[n[1]][j][q] * L[n[2]][k][r] * w[p][q][r]
                for p in range(l[0]):
                    for q in range(l[1]):
                        for r in range(l[2]):
                            R[i][j][k][p][q][r] = (
                                N[n[0]][i][p] * M[n[1]][j][q] * L[n[2]][k][r] * w[p][q][r]) / a[i][j][k]
    return R


def trans(a, b, c):
    d = np.array([[a],
                  [b],
                  [c]])
    return d


def rotx(a):
    b = np.array([[1, 0,            0],
                  [0, math.cos(a),  - math.sin(a)],
                  [0, math.sin(a),  math.cos(a)]])
    return b


def roty(a):
    b = np.array([[math.cos(a),     0, math.sin(a)],
                  [0,               1, 0],
                  [- math.sin(a),   0, math.cos(a)]])
    return b


def rotz(a):
    b = np.array([[math.cos(a), - math.sin(a),  0],
                  [math.sin(a), math.cos(a),    0],
                  [0,           0,              1]])
    return b


def stretchxyz(a, b, c):
    d = np.array([[a, 0, 0],
                  [0, b, 0],
                  [0, 0, c]])
    return d


def shearx(a):
    b = np.array([[1, math.tan(a),  0],
                  [0, 1,            0],
                  [0, 0,            1]])
    return b


def sheary(a):
    b = np.array([[1,           0, 0],
                  [math.tan(a), 1, 0],
                  [0,           0, 1]])
    return b


def affine_transformation_2D(CP, l, stretch_x, stretch_y, stretch_z, trans_x, trans_y, trans_z, theta_x, theta_y, theta_z, shear_x, shear_y):
    CP_3D = np.zeros((l, 3))
    for i in range(l):
        for j in range(3):
            if j == 2:
                CP_3D[i][j] = 0
            else:
                CP_3D[i][j] = CP[i][j]
    CP_new = np.zeros((l, 3))
    rotxyz = np.dot(np.dot(rotx(theta_x), roty(theta_y)), rotz(theta_z))
    shearxy = np.dot(shearx(shear_x), sheary(shear_y))

    for i in range(l):
        a = np.dot(np.dot(np.dot(shearxy, rotxyz), stretchxyz(
            stretch_x, stretch_y, stretch_z)), (CP_3D[i, :]).reshape(-1, 1))
        b = trans(trans_x, trans_y, trans_z)
        CP_new[i][0] = (a + b)[0][0]
        CP_new[i][1] = (a + b)[1][0]
        CP_new[i][2] = (a + b)[2][0]
    CP_3D = CP_new
    for i in range(l):
        for j in range(2):
            CP[i][j] = CP_3D[i][j]
    return CP


def affine_transformation_3D(CP, l, stretch_x, stretch_y, stretch_z, trans_x, trans_y, trans_z, theta_x, theta_y, theta_z, shear_x, shear_y):
    CP_3D = np.zeros((l, 3))
    for i in range(l):
        for j in range(3):
            CP_3D[i][j] = CP[i][j]
    CP_new = np.zeros((l, 3))
    rotxyz = np.dot(np.dot(rotx(theta_x), roty(theta_y)), rotz(theta_z))
    shearxy = np.dot(shearx(shear_x), sheary(shear_y))

    for i in range(l):
        a = np.dot(np.dot(np.dot(shearxy, rotxyz), stretchxyz(
            stretch_x, stretch_y, stretch_z)), (CP_3D[i, :]).reshape(-1, 1))
        b = trans(trans_x, trans_y, trans_z)
        CP_new[i][0] = (a + b)[0][0]
        CP_new[i][1] = (a + b)[1][0]
        CP_new[i][2] = (a + b)[2][0]
    CP = CP_new
    return CP


def get_4element_2D(Sx_vec , Sy_vec, Sz_vec, point_vec, a, b):
    for i in range(a):
        for j in range(b):
            point_vec[i,j,0,:] = np.array([Sx_vec[i+1,j], Sy_vec[i+1,j], Sz_vec[i+1,j]])
            point_vec[i,j,1,:] = np.array([Sx_vec[i,j], Sy_vec[i,j], Sz_vec[i,j]])
            point_vec[i,j,2,:] = np.array([Sx_vec[i,j+1], Sy_vec[i,j+1], Sz_vec[i,j+1]])
            point_vec[i,j,3,:] = np.array([Sx_vec[i+1,j+1], Sy_vec[i+1,j+1], Sz_vec[i+1,j+1]])
    return point_vec


def get_4element_3D(Sx_vec, Sy_vec, Sz_vec, point_vec_i1, point_vec_i2, point_vec_j1, point_vec_j2, point_vec_k1, point_vec_k2, a, b, c):
    point_vec_i1 = get_4element_2D(Sx_vec[:,:,0], Sy_vec[:,:,0], Sz_vec[:,:,0], point_vec_i1, a, b)
    point_vec_i2 = get_4element_2D(Sx_vec[:,:,-1], Sy_vec[:,:,-1], Sz_vec[:,:,-1], point_vec_i2, a, b)
    point_vec_j1 = get_4element_2D(Sx_vec[0,:,:], Sy_vec[0,:,:], Sz_vec[0,:,:], point_vec_j1, b, c)
    point_vec_j2 = get_4element_2D(Sx_vec[-1,:,:], Sy_vec[-1,:,:], Sz_vec[-1,:,:], point_vec_j2, b, c)
    point_vec_k1 = get_4element_2D(Sx_vec[:,0,:], Sy_vec[:,0,:], Sz_vec[:,0,:], point_vec_k1, a, c)
    point_vec_k2 = get_4element_2D(Sx_vec[:,-1,:], Sy_vec[:,-1,:], Sz_vec[:,-1,:], point_vec_k2, a, c)
    return point_vec_i1, point_vec_i2, point_vec_j1, point_vec_j2, point_vec_k1, point_vec_k2


def transform_triangle_2D(a, b, point_vec, mesh, n_vec):
    for i in range(a):
        for j in range(b):
            for k in range(3):
                mesh[i][j][0][0][k] = point_vec[i][j][0][k]
                mesh[i][j][0][1][k] = point_vec[i][j][1][k]
                mesh[i][j][0][2][k] = point_vec[i][j][2][k]
            for k in range(3):
                mesh[i][j][1][0][k] = point_vec[i][j][2][k]
                mesh[i][j][1][1][k] = point_vec[i][j][3][k]
                mesh[i][j][1][2][k] = point_vec[i][j][0][k]
    p = np.zeros((2, 3))
    for i in range(a):
        for j in range(b):
            for k in range(2):
                for l in range(3):
                    p[0][l] = mesh[i][j][k][0][l] - mesh[i][j][k][1][l]
                    p[1][l] = mesh[i][j][k][1][l] - mesh[i][j][k][2][l]
                    q = np.cross(p[0,:], p[1,:])
                    r = np.linalg.norm(q)
                    if r != 0:
                        n_vec[i,j,k,:] = q / r
                    if r == 0:
                        n_vec[i,j,k,0] = 1.0
                        n_vec[i,j,k,1] = 0.0
                        n_vec[i,j,k,2] = 0.0
    return mesh, n_vec


def write_stl_header(solid_name):
    solid_name_stl = solid_name + ".stl"
    f = open(solid_name_stl, 'w')
    f.write('solid ')
    f.write(solid_name)
    f.write('\n')


def write_stl_main(solid_name, a, b, mesh, n_vec):
    solid_name_stl = solid_name + ".stl"
    f = open(solid_name_stl, 'a')
    for i in range(a):
        for j in range(b):
            for k in range(2):
                f.write('   facet nomal ')
                f.write(str(n_vec[i,j,k,0]))
                f.write(' ')
                f.write(str(n_vec[i,j,k,1]))
                f.write(' ')
                f.write(str(n_vec[i,j,k,2]))
                f.write('\n')
                f.write('       outer loop\n')
                f.write('           vertex ')
                f.write(str(mesh[i][j][k][0][0]))
                f.write(' ')
                f.write(str(mesh[i][j][k][0][1]))
                f.write(' ')
                f.write(str(mesh[i][j][k][0][2]))
                f.write('\n')
                f.write('           vertex ')
                f.write(str(mesh[i][j][k][1][0]))
                f.write(' ')
                f.write(str(mesh[i][j][k][1][1]))
                f.write(' ')
                f.write(str(mesh[i][j][k][1][2]))
                f.write('\n')
                f.write('           vertex ')
                f.write(str(mesh[i][j][k][2][0]))
                f.write(' ')
                f.write(str(mesh[i][j][k][2][1]))
                f.write(' ')
                f.write(str(mesh[i][j][k][2][2]))
                f.write('\n')
                f.write('      endlooop\n')
                f.write('   endfacet\n')


def write_stl_footer(solid_name):
    solid_name_stl = solid_name + ".stl"
    f = open(solid_name_stl, 'a')
    f.write('endsolid')
    f.close()


def make_stl_2D(solid_name, delta, Sx_vec, Sy_vec, Sz_vec):
    a = int(delta[0] - 1)
    b = int(delta[1] - 1)
    point_vec = np.zeros((a, b, 4, 3))
    point_vec = get_4element_2D(Sx_vec, Sy_vec, Sz_vec, point_vec, a, b)
    mesh = np.zeros((a, b, 2, 3, 3))
    n_vec = np.zeros((a, b, 2, 3))
    mesh, n_vec = transform_triangle_2D(a, b, point_vec, mesh, n_vec)
    write_stl_header(solid_name)
    write_stl_main(solid_name, a, b, mesh, n_vec)
    write_stl_footer(solid_name)


def make_stl_3D(solid_name, delta, Sx_vec, Sy_vec, Sz_vec):
    a = int(delta[0] - 1)
    b = int(delta[1] - 1)
    c = int(delta[2] - 1)
    point_vec_i1 = np.zeros((a, b, 4, 3))
    point_vec_i2 = np.zeros((a, b, 4, 3))
    point_vec_j1 = np.zeros((b, c, 4, 3))
    point_vec_j2 = np.zeros((b, c, 4, 3))
    point_vec_k1 = np.zeros((a, c, 4, 3))
    point_vec_k2 = np.zeros((a, c, 4, 3))
    point_vec_i1, point_vec_i2, point_vec_j1, point_vec_j2, point_vec_k1, point_vec_k2 = get_4element_3D(Sx_vec, Sy_vec, Sz_vec, point_vec_i1, point_vec_i2, point_vec_j1, point_vec_j2, point_vec_k1, point_vec_k2, a, b, c)
    mesh_i1 = np.zeros((a, b, 2, 3, 3))
    n_vec_i1 = np.zeros((a, b, 2, 3))
    mesh_i2 = np.zeros((a, b, 2, 3, 3))
    n_vec_i2 = np.zeros((a, b, 2, 3))
    mesh_j1 = np.zeros((b, c, 2, 3, 3))
    n_vec_j1 = np.zeros((b, c, 2, 3))
    mesh_j2 = np.zeros((b, c, 2, 3, 3))
    n_vec_j2 = np.zeros((b, c, 2, 3))
    mesh_k1 = np.zeros((a, c, 2, 3, 3))
    n_vec_k1 = np.zeros((a, c, 2, 3))
    mesh_k2 = np.zeros((a, c, 2, 3, 3))
    n_vec_k2 = np.zeros((a, c, 2, 3))
    mesh_i1, n_vec_i1 = transform_triangle_2D(a, b, point_vec_i1, mesh_i1, n_vec_i1)
    mesh_i2, n_vec_i2 = transform_triangle_2D(a, b, point_vec_i2, mesh_i2, n_vec_i2)
    mesh_j1, n_vec_j1 = transform_triangle_2D(b, c, point_vec_j1, mesh_j1, n_vec_j1)
    mesh_j2, n_vec_j2 = transform_triangle_2D(b, c, point_vec_j2, mesh_j2, n_vec_j2)
    mesh_k1, n_vec_k1 = transform_triangle_2D(a, c, point_vec_k1, mesh_k1, n_vec_k1)
    mesh_k2, n_vec_k2 = transform_triangle_2D(a, c, point_vec_k2, mesh_k2, n_vec_k2)
    write_stl_header(solid_name)
    if np.any(mesh_i1 != mesh_i2):
        write_stl_main(solid_name, a, b, mesh_i1, n_vec_i1)
        write_stl_main(solid_name, a, b, mesh_i2, n_vec_i2)
    if np.any(mesh_j1 != mesh_j2):
        write_stl_main(solid_name, b, c, mesh_j1, n_vec_j1)
        write_stl_main(solid_name, b, c, mesh_j2, n_vec_j2)
    if np.any(mesh_k1 != mesh_k2):
        write_stl_main(solid_name, a, c, mesh_k1, n_vec_k1)
        write_stl_main(solid_name, a, c, mesh_k2, n_vec_k2)
    write_stl_footer(solid_name)