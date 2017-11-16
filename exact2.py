#!/usr/bin/env python
import numpy as np
import math
import time

def cispace(d,occ):
    #d: dim of spin-orbital
    #occ: num of electrons
    #return all bases of CI (in ONV)

    if d < occ or occ < 0 or d <= 0:
        return []
    if d == 1:
            return [[occ]]
    else:
        ans1 = [i + [1] for i in cispace(d-1,occ-1)]
        ans2 = [i + [0] for i in cispace(d-1,occ)]
        return ans1 + ans2

def xorloc(a,b):
    #a, b: 2 bases of CI
    #return the location where a is 1 while b is 0

    d = len(a)
    ans = []
    for i in range(d):
        if a[i] == 1 and b[i] == 0:
            ans = ans + [i]
    return ans

def xordist(a,b):
    #a, b: lists
    #return Hamming distance

    l1 = np.array(a)
    l2 = np.array(b)
    return np.shape(np.nonzero(l1-l2)[0])[0]

def e1s(m,n):
    #in ONV, odd: alpha spin; even: beta spin
    #return <m|h|n>

    e1 = np.load("h1e.npy")
    ans = 0.0
    if m%2 == n%2:
        ans += e1[m//2][n//2]
    return ans

def e2s(m,n,p,q):
    #in ONV, odd: alpha spin; even: beta spin
    #return <mn||pq>
    
    e2 = np.load("h2e.npy")
    ans = 0.0
    if m%2 == p%2 and n%2 == q%2:
        ans += e2[m//2][p//2][n//2][q//2]
    if m%2 == q%2 and n%2 == p%2:
        ans -= e2[m//2][q//2][n//2][p//2]
    return ans

def phase(a, ldiff):
    #a: basis of CI
    #ldiff: whichi list to put in the very beginning

    ans = 0
    for i in range(len(ldiff)):
        for j in range(ldiff[i]):
            if a[j] == 1:
                ans += 1
        ans -= i
        ans = ans % 2
    return (-1)**ans

def slater(a,b):
    #a, b: 2 bases of CI (in ONV, odd: alpha spin; even: beta spin)

    dist = xordist(a,b)
    d = len(a)
    if dist == 0:
        ans = 0.0
        for i in range(d):
            if a[i] == 1:
                ans += e1s(i,i)
                for j in range(i):
                    if a[j] == 1:
                        ans += e2s(i,j,i,j)
        return ans

    elif dist == 2:
        m = xorloc(a,b)
        n = xorloc(b,a)
        ans = e1s(m[0],n[0])
        for i in range(d):
            if a[i] == 1 and i != m[0] and i != n[0]:
                ans += e2s(m[0],i,n[0],i)
        return ans*phase(a,m)*phase(b,n)

    elif dist == 4:
        m = xorloc(a,b)
        n = xorloc(b,a)
        return e2s(m[0],m[1],n[0],n[1])*phase(a,m)*phase(b,n)

    else:
        return 0.0

t1 = time.time()

e1 = np.load("h1e.npy")
d = e1.shape[0]
s = 0
index = []
alpha = cispace(d, (d+s)/2)
beta = cispace(d, (d-s)/2)

for a in alpha:
    for b in beta:
        tmp = []
        for i in range(d):
            tmp += [a[i]]
            tmp += [b[i]]
        index += [tmp]

dim = len(index)
ham = np.zeros((dim,dim))

for i in range(dim):
    print i
    for j in range(i+1):
        ham[i][j] = slater(index[i],index[j])
        ham[j][i] = ham[i][j]

en = min(np.linalg.eig(ham)[0])

t2 = time.time()

print en
print t2-t1

