#!/usr/bin/env python
import numpy as np
import math
import time

def foo(d,occ):
    #d: dim of spin-orbital
    #occ: num of electrons
    #return all bases of CI (in ONV, odd: alpha spin; even: beta spin)

    if d < occ or occ < 0 or d <= 0:
        return []
    if d == 1:
            return [[occ]]
    else:
        ans1 = [i + [1] for i in foo(d-1,occ-1)]
        ans2 = [i + [0] for i in foo(d-1,occ)]
        return ans1 + ans2

def xorloc1(a,b):
    #a, b: 2 bases of CI
    #return the location where a is 1 while b is 0

    d = len(a)
    ans = []
    for i in range(d):
        if a[i] == 1 and b[i] == 0:
            ans = ans + [i]
    return ans

def xorloc2(a,b):
    #a, b: 2 bases of CI
    #return the location where b is 1 while a is 0

    d = len(a)
    ans = []
    for i in range(d):
        if b[i] == 1 and a[i] == 0:
            ans = ans + [i]
    return ans

def xordist(a,b):
    #a, b: lists
    #return Hamming distance

    l1 = np.array(a)
    l2 = np.array(b)
    return np.shape(np.nonzero(l1-l2)[0])[0]

def e1s(m,n):
    #return <m|h|n>

    e1 = np.load("h1e.npy")
    ans = 0.0
    if m%2 == n%2:
        ans += e1[m//2][n//2]
    return ans

def e2s(m,n,p,q):
    #return <mn||pq>
    
    e2 = np.load("h2e.npy")
    ans = 0.0
    if m%2 == p%2 and n%2 == q%2:
        ans += e2[m//2][p//2][n//2][q//2]
    if m%2 == q%2 and n%2 == p%2:
        ans -= e2[m//2][q//2][n//2][p//2]
    return ans

def phase1(a,b):
    #a, b: 2 bases of CI

    da = xorloc1(a,b)[0]
    db = xorloc2(a,b)[0]
    ans = 0
    if da < db:
        a, b = b, a
        da, db = db, da
    for i in range(db+1,da):
        if a[i] == 1:
            ans += 1
    return (-1)**ans

def phase2(a,b):
    #a, b: 2 bases of CI

    da = xorloc1(a,b)[0]
    db = xorloc2(a,b)[0]
    ans = 0
    if da < db:
        a, b = b, a
        da, db = db, da
    for i in range(db+1,da):
        if a[i] == 1:
            ans += 1
    da = xorloc1(a,b)[1]
    db = xorloc2(a,b)[1]
    if da > db:
        a, b = b, a
        da, db = db, da
    for i in range(da+1,db):
        if a[i] == 1:
            ans += 1
    return (-1)**ans

def slater(a,b):
    #a, b: 2 bases of CI

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
        m = xorloc1(a,b)[0]
        n = xorloc2(a,b)[0]
        ans = e1s(m,n)
        for i in range(d):
            if a[i] == 1 and i != m and i != n:
                ans += e2s(m,i,n,i)
        return ans*phase1(a,b)

    elif dist == 4:
        m = xorloc1(a,b)[0]
        n = xorloc1(a,b)[1]
        p = xorloc2(a,b)[0]
        q = xorloc2(a,b)[1]
        return e2s(m,n,p,q)*phase2(a,b)

    else:
        return 0.0
