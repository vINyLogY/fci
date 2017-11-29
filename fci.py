import numpy as np
import math

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

def merge(alpha,beta):
    index = []
    for a in alpha:
        for b in beta:
            tmp = []
            for i in range(len(a)):
                tmp += [a[i]]
                tmp += [b[i]]
            index += [tmp]
    return index

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
        #rank = np.linalg.matrix_rank(V.T.dot(V))
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
        daVec = resVec / (ritzVal - Da + 1.0e-5)    #avoid divided by zero
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

    return ritzVal
