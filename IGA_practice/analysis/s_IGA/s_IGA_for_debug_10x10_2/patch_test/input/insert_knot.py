import numpy as np
import math

# PI = math.pi
# x1 = math.tan(PI / 8.)
# x2 = math.sin(PI / 4.)
# x3 = math.cos(PI / 8.)

# print("{:.20e}".format(x1))
# print("{:.20e}".format(x2))
# print("{:.20e}".format(x3))

print("---------------")

a = 2 # 挿入数

A = np.zeros((a))

for i in range(a):
    A[i] = (i + 1.) / (a + 1.)
    print("{:.20e}".format(A[i]))