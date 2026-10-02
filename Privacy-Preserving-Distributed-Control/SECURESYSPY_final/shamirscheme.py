import numpy as np
#np.random.seed(1)
# Creates shares of secrets using Shamir's secret sharing scheme.(input can be a float)
#using share(F, x, t, n) function a large finite field must be defined (1000 times larger than biggest used number in the field)
#finite field considering negative nmbers
def truncate(num, n):
    integer = int(num * (10**n))/(10**n)
    return float(integer)

#return the float number in such a way that each element of finfield related to  0,....(p-1)/2,-(p-1)/2,....,-1
def F2R(F,a,scale):
    p=int(F.cardinality())
    if int(a) <=  (p-1)/2:
        return float(int(a)/10**scale)
    elif int(a)>(p-1)/2:
        return float(-(p-int(a))/10**scale)
# return corresponidng integer in field F with scaling
def R2F(F,a,scale):
    a=float(a)
    a=truncate(a, scale)
    a=F(int(a*10**scale))   
    return F(a)

def share(F, x, t, n,scale):
    x=R2F(F,x,scale)
    shares = []
    c = [F.random_element() for i in range(t)]
    for i in range(1, n + 1):
        s = x
        for j in range(1, t + 1):
            s += c[j - 1] * F(i)**j
        shares.append(s)
    s = np.array(shares)
    return s

#used in mult
def intshare(F, x, t, n):
    x=F(x)
    shares = []
    c = [F.random_element() for i in range(t)]
    for i in range(1, n + 1):
        s = x
        for j in range(1, t + 1):
            s += c[j - 1] * F(i)**j
        shares.append(s)
    s = np.array(shares)
    return s
# reconstruct secret.#r is recombination vector
# Creates the "recombination"-vector used to reconstruct a secret from its shares.
def basispoly(F,n):
    r = []
    C = range(1, n + 1)
    for i in range(1, n + 1):
        c = [k for k in C if k != i]
        p = 1
        for j in range(n - 1):
            p *= -F(c[j]) / (F(i) - F(c[j]))
        r.append(p)
    return r

def rec(F,x,r,scale):
    res = F(0)
    n = len(x)
    for i in range(len(x)):
        res += x[i] * r[i]
    return F2R(F,res,scale)

def intrec(F,x,r):
    res = F(0)
    n = len(x)
    for i in range(len(x)):
        res += x[i] * r[i]
    return res