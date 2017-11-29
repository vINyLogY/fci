#!/usr/bin/env python
import numpy as np
import math
import time
from fci1 import *

#input
e1 = np.load("h1e.npy")
e2 = np.load("h2e.npy")
d = e1.shape[0]
ne = d  #ne = alpha + beta
s = 0   #s = alpha - beta

#form cispace string
alpha = cispace(d, (ne+s)/2)
beta = cispace(d, (ne-s)/2)
dim = len(alpha)*len(beta)
space = merge(alpha, beta)

#form DA
Da = np.zeros(dim)
for i in range(dim):
    C = np.zeros(dim)
    C[i] = 1.0
    Da[i] = C.dot(direct(C,space,e1,e2))

#form initial search space
t1 = time.time()
V1 = np.zeros(dim)
V1[0] = 1.0

#davidson algorithm
HProd = lambda x: direct(x,space,e1,e2)
ans = davidson(HProd,Da,V1,1.0e-4,10)
t2 = time.time()

#output
t = t2 - t1
print("ans = %.8f" % ans)
print("time = %f" % t)
