import numpy as np
from crc import count_crc

def polar_encode(u, F, N):
    c = np.zeros(N)
    c[F==0] = u
    x = polar_encode_recursive(c.copy())
    return x, c

def polar_encode_crc(u, F, N, poly):
    c = np.zeros(N)
    T = count_crc(u, poly)
    u_T = np.concatenate((u, T))
    c[F==0] = u_T
    x = polar_encode_recursive(c.copy())
    return x, c

def polar_encode_recursive(c):
    N = len(c)
    c = c.astype(int)
    if N == 2:
        return np.array([c[0]^c[1], c[1]])
    x1 = polar_encode_recursive(c[:int(N/2)])
    x2 = polar_encode_recursive(c[int(N/2):])
    return np.concatenate((x1^x2, x2))