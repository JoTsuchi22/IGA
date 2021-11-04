import numpy as np

a = 11 # 挿入数

A = np.zeros((a))

for i in range(a):
    A[i] = (i + 1.) / (a + 1.)
    print("{:.20e}".format(A[i]))