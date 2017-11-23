#!/usr/bin/env python
import numpy as np
import math
import time
from numpy.random import rand, randint
from scipy.sparse.linalg import eigsh
from scipy.sparse import csc_matrix

def minIndex(li):
    ans = 0
    for i in range(len(li)):
        if li[i] < li[ans]:
            ans = i
    return ans        

def davidson(H,V,err,maxIter):
    #initialization
    dim = H.shape[0]
    numIter = 0
    last = 0.0
    Da = H.diagonal()
    HV = H.dot(V)
    A = V.T.dot(HV)

    #iteration
    while numIter <= maxIter:
        #print(20*"*")
        if numIter != 0:
            last = ritzVal
        numIter += 1
        #print("Iter: %d" % numIter)
        rank = np.linalg.matrix_rank(V.T.dot(V))
        #print("Rank of V*V: %d" % rank)

        
        #form ritz value and vector
        eigensolver = np.linalg.eig(A)
        ritzIndex = minIndex(eigensolver[0])
        ritzVal = eigensolver[0][ritzIndex]
        rVec = eigensolver[1][ritzIndex]
        #print("Ritz Value: %f" % (ritzVal))

        #check convergence
        resVec = HV.dot(rVec) - ritzVal * V.dot(rVec)
        conv = ritzVal - last
        conv1 = np.linalg.norm(resVec)
        #print("norm: %f" % conv1)
        if conv < err**2 and conv > 0.0 - err**2:
            break
        if conv1 < err:
            break

        #expand the search space
        daVec = resVec / (ritzVal - Da + 1.0e-5)
        for i in V.T:
            daVec -= i.dot(daVec) * i
        norm = np.linalg.norm(daVec)
        if norm < err:
            break
        #print("norm of daVec %f" % norm)
        daVec /= norm
        Vi = daVec.reshape(dim,1)
        HVi = H.dot(Vi)
        Ai = V.T.dot(HVi)
        ai = Vi.T.dot(HVi)
        A = np.vstack([np.hstack([A, Ai]), np.hstack([Ai.T, ai])]) 
        V = np.hstack([V, Vi])
        HV = np.hstack([HV, HVi])

    return ritzVal, conv1, numIter, V, rank

def randmat(dim):
    A = np.zeros((dim,dim))
    for i in range(dim):
        for j in range(i+1):
            if rand() < 1.0 / dim:
                A[i][j] = rand()
                A[j][i] = A[i][j]
    return A

frac = 0.5
for i in range(3000):
    A = randmat(3000)
    A = csc_matrix(A)
    t0 = time.time()
    ep = eigsh(A)
    ans_ref = ep[0][0]
    t1 = time.time()
    V1 = math.sqrt(frac) * (ep[1].T)[0:1].T + math.sqrt(1.0-frac) * (ep[1].T)[1:2].T
    ans, norm, numIter, V, rank = davidson(A,V1,1.0e-6,100)
    t2 = time.time()
    t = t2 - t1
    t_ref = t1 - t0
    err = (ans-ans_ref)
    rt = t/t_ref
    print("Error = %f" % err)
    #print("relative time = %f" % rt)
    if err > 1.0 or err < -1.0:
        print(20*"-")
        print(ans)
        print("time = %f" % t)
        print("norm = %f" % norm)
        print("num of iter = %d" % numIter)
        print("Rank of V*V: %d" % rank)
        print("reference ans = %f" % ans_ref)
        print("reference time = %f" % t_ref)
        print(20*"-")
        break
