import numpy as np

def count_crc(u, poly):
    N = len(poly)
    T = np.zeros(N)
    u_0 = np.concatenate((u, np.zeros(N)))
    for i in range(len(u_0)):
        T0 = T[0]
        T = np.concatenate((T[1:], [u_0[i]])).astype(int)
        if T0 == 1:
            T = T ^ poly
    return T

def check_crc(u, poly):
    N = len(poly)
    T = np.zeros(N)
    for i in range(len(u)):
        T0 = T[0]
        T = np.concatenate((T[1:], [u[i]])).astype(int)
        if T0 == 1:
            T = T ^ poly
    return T.sum()