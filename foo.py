#!/usr/bin/env python
import numpy as np
import math
import time

def minIndex(li):
    ans = 0
    for i in range(len(li)):
        if li[i] < li[ans]:
            ans = i
    return ans        

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
            ans += [i]
    return ans

def xordist(a,b):
    #a, b: lists
    #return Hamming distance

    l1 = np.array(a)
    l2 = np.array(b)
    return np.shape(np.nonzero(l1-l2)[0])[0]

def e1s(m,n,e1):
    #in ONV, odd: alpha spin; even: beta spin
    #return <m|h|n>

    ans = 0.0
    if m%2 == n%2:
        ans += e1[m//2][n//2]
    return ans

def e2s(m,n,p,q,e2):
    #in ONV, odd: alpha spin; even: beta spin
    #return <mn||pq>
    
    ans = 0.0
    if m%2 == p%2 and n%2 == q%2:
        ans += e2[m//2][p//2][n//2][q//2]
    if m%2 == q%2 and n%2 == p%2:
        ans -= e2[m//2][q//2][n//2][p//2]
    return ans

def phase(a,ldiff):
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

def slater(a,b,e1,e2):
    #a, b: 2 bases of CI (in ONV, odd: alpha spin; even: beta spin)

    dist = xordist(a,b)
    d = len(a)
    if dist == 0:
        ans = 0.0
        for i in range(d):
            if a[i] == 1:
                ans += e1s(i,i,e1)
                for j in range(i):
                    if a[j] == 1:
                        ans += e2s(i,j,i,j,e2)
        return ans

    elif dist == 2:
        m = xorloc(a,b)
        n = xorloc(b,a)
        ans = e1s(m[0],n[0],e1)
        for i in range(d):
            if a[i] == 1 and i != m[0] and i != n[0]:
                ans += e2s(m[0],i,n[0],i,e2)
        return ans*phase(a,m)*phase(b,n)

    elif dist == 4:
        m = xorloc(a,b)
        n = xorloc(b,a)
        return e2s(m[0],m[1],n[0],n[1],e2)*phase(a,m)*phase(b,n)

    else:
        return 0.0

def davidson(H,V,err,maxIter):
    #initialization
    dim =len(H)
    numIter = 0
    last = 0.0
    Da = np.diag(H)
    HV = H.dot(V)
    A = V.T.dot(HV)
    Orth = np.eye(dim)
    for i in V.T:
        Vi = i.reshape(dim,1)
        Orth = Orth.dot(np.eye(dim)-Vi.dot(Vi.T))

    #iteration
    while numIter <= maxIter:
        #print("----------")
        #print("Iter: %d" % numIter)
        if numIter != 0:
            last = ritzVal
        numIter += 1

        #form ritz value and ritz vector
        eigensolver = np.linalg.eig(A)
        ritzIndex = minIndex(eigensolver[0])
        ritzVal = eigensolver[0][ritzIndex]
        rVec = eigensolver[1][ritzIndex]
        #print("Ritz Value: %f" % (ritzVal))

        #check convergence
        resVec = HV.dot(rVec) - ritzVal * V.dot(rVec)
        conv = ritzVal - last
        conv1 = np.linalg.norm(resVec)
        if conv < err and conv > 0.0 - err:
            break
        #print("conv: %f" % conv)

        #expand the search space
        daVec = resVec / (ritzVal - Da)
        daVec = Orth.dot(daVec)
        daVec /= np.linalg.norm(daVec)
        Vi = daVec.reshape(dim,1)
        HVi = H.dot(Vi)
        Ai = V.T.dot(HVi)
        ai = Vi.T.dot(HVi)
        Orth = Orth.dot(np.eye(dim)-Vi.dot(Vi.T))
        A = np.vstack([np.hstack([A, Ai]), np.hstack([Ai.T, ai])]) 
        V = np.hstack([V, Vi])
        HV = np.hstack([HV, HVi])

    return ritzVal

#input
e1 = np.load("h1e.npy")
e2 = np.load("h2e.npy")
d = e1.shape[0]
ne = d  #ne = alpha + beta
s = 0   #s = alpha - beta

# form cispace string
index = []
alpha = cispace(d, (ne+s)/2)
beta = cispace(d, (ne-s)/2)
for a in alpha:
    for b in beta:
        tmp = []
        for i in range(d):
            tmp += [a[i]]
            tmp += [b[i]]
        index += [tmp]
dim = len(index)

#form Hamitonian
Ham = np.zeros((dim,dim))
for i in range(dim):
    for j in range(i+1):
        Ham[i][j] = slater(index[i],index[j],e1,e2)
        Ham[j][i] = Ham[i][j]
en = min(np.linalg.eig(Ham)[0])

#form initial search space
salpha = cispace((ne+s)/2+1, (ne+s)/2)
salpha = [i + (d - len(i)) * [0] for i in salpha]
sbeta = cispace((ne-s)/2+1, (ne-s)/2)
sbeta = [i + (d - len(i)) * [0] for i in sbeta]
sindex = []
for a in salpha:
    for b in sbeta:
        tmp = []
        for i in range(d):
            tmp += [a[i]]
            tmp += [b[i]]
        sindex += [tmp]
sdim = len(sindex)
V = np.zeros((dim,sdim))
for i in range(sdim):
    V[index.index(sindex[i])][i] = 1.0
V1 = np.zeros((dim,3))
for i in range(3):
    V1[-i-1][i] = 1.0

#davidson algorithm
ans = davidson(Ham, V1, 1.0e-9, dim)

#output
print en
print ans

