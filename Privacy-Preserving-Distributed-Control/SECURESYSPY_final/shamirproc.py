#shamirproc
import numpy as np
import shamirscheme as sss
import math

def add(F,a,b,t,n,scale):
    r =sss.basispoly(F,n) 
    if 'numpy' in str(type(a)):
        sa=a
        sb=b
    else:
        sa=sss.share(F,a,t,n,scale)
        sb=sss.share(F,b,t,n,scale)
    Addition=sa+sb
    add.Addition=Addition
    return np.array(Addition)

def Idiv(F,q,k,n,t,r):#q/k
    #party 1
    rnd=F.random_element()
    rndonk=math.floor(int(rnd)/k)
    rndshares=sss.intshare(F, rnd, t, n)
    rndonkshares=sss.intshare(F, rndonk, t, n)
    #parties compute shares of Z
    z=q+rndshares
    #open z for P2
    z=sss.intrec(F,z,r)
    # P2 computes
    zonk=math.floor(int(z)/k)
    zonkshares = sss.intshare(F, zonk, t, n)
    y = zonkshares-rndonkshares
    return y

#mult a,b which have been shared with scaling
def mult(F,a,b,t,n,scale):
    r =sss.basispoly(F,n)    
    if 'numpy' in str(type(b)) and str(type(a)):
        sa=Idiv(F,a,10**(scale/2),n,t,r)
        sb=Idiv(F,b,10**(scale/2),n,t,r)
        sh=sa*sb
    else:
        sa=sss.share(F,a,t,n,scale/2)
        sb=sss.share(F,b,t,n,scale/2)
        sh=sa*sb
    sasbshares = []
    #abshares[i,;] is the row correspond to shares of ki=sa[i]*sb[i] which are hi(x) i=1,...,n
    for i in range(n):
        sasbshares.append((sss.intshare(F,sh[i],t,n)))
    abshares = []
    for i in range(n): 
        yj=F(0)
        for j in range(n):
            yj += sasbshares[j][i]* r[j]
        abshares.append(yj)
    mult.multiplication=abshares
    return np.array(abshares)

#F=GF(next_prime(2^100))
#a=-5.256
#b=321
#t=1
#n=3
#print(mult(F,a,b,t,n))
#print(add(F,a,b,t,n))