#!/usr/bin/env python
import numpy as np
import math
import time
from fci import *

#input
e1 = np.load("h1e.npy")
e2 = np.load("h2e.npy")
d = e1.shape[0]
ne = d  #ne = alpha + beta
s = 0   #s = alpha - beta

#form cispace string
alpha = cispace(d, (ne+s)/2)
beta = cispace(d, (ne-s)/2)
index = merge(alpha, beta)
dim = len(index)

#form Hamitonian
Ham = np.zeros((dim,dim))
for i in range(dim):
    for j in range(i+1):
        Ham[i][j] = slater(index[i],index[j],e1,e2)
        Ham[j][i] = Ham[i][j]

#form initial search space
t1 = time.time()
V1 = np.zeros((dim,1))
for i in range(1):
    V1[-i-1][i] = 1.0

#davidson algorithm
ans = davidson(Ham, V1, 1.0e-6, dim)
t2 = time.time()

#output
en = min(np.linalg.eig(Ham)[0])
t3 = time.time()
t = t2 - t1
tr = t3 - t2
print("ans = %.8f" % ans)
print("time = %f" % t)
print("reference ans = %.8f" % en)
print("reference time = %f" % tr)

