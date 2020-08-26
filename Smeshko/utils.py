import numpy as np


def f(a,b):
    return np.sign(a*b)*min(np.abs(a), np.abs(b))

def f1(a,b):
    return np.sign(a*b)*np.min(np.abs(np.stack([a, b])), axis=0)

def g(a,b,c):
    if c == 0:
        return a + b
    else:
        return b - a
    
def g1(a,b,c):
    return b  + (1 - 2 * c) * a

def f_N(llr):
    N = len(llr)
    llr_1 = llr[:int(N/2)]
    llr_2 = llr[int(N/2):]
    return [f(a,b) for a,b in zip(llr_1, llr_2)] # f(llr_1, llr_2) #

def f_N1(llr):
    N = len(llr)
    llr_1 = llr[:int(N/2)]
    llr_2 = llr[int(N/2):]
    return f1(llr_1, llr_2) #


def g_N(llr, u):
    N = len(llr)
    llr_1 = llr[:int(N/2)]
    llr_2 = llr[int(N/2):]
    return [g(a,b,c) for a,b,c in zip(llr_1, llr_2, u)]

def g_N1(llr, u):
    N = len(llr)
    llr_1 = llr[:int(N/2)]
    llr_2 = llr[int(N/2):]
    return g1(llr_1, llr_2, u)

def check(L, F):
    if F == 1:
        return 0
    if L >= 0:
        return 0
    else:
        return 1
    
def check1(L, F):
    if F == 1 or L >= 0:
        return 0
    return 1

def check_list(Li, F):
    if F == 1:
        if Li >= 0:
            return [[0, 0]]
        return [[0, np.abs(Li)]]
    if Li >= 0:
        return [[0, 0], [1, np.abs(Li)]]
    return [[1, 0], [0, np.abs(Li)]]
