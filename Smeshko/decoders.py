import numpy as np
from utils import check_list, f, f_N, g, g_N, f_N1, g_N1, f1, g1, check, check1
from structures import StructureForList, StructureForStack
from crc import check_crc

def polar_decode_sc1(in_llr, F, u):
    n = len(in_llr) // 2
    if n == 1:
        L_1 = f(in_llr[0], in_llr[1])
        u_1 = check1(L_1, F[0])
        L_2 = g(in_llr[0], in_llr[1], u_1)
        u_2 = check1(L_2, F[1])
        u.append(u_1)
        u.append(u_2)
        return u, in_llr, np.array([u_1^u_2, u_2])
    else:
        u, out_llr1, x1 = polar_decode_sc1(f_N1(in_llr), F[:n], u)
        u, out_llr2, x2 = polar_decode_sc1(g_N1(in_llr, x1), F[n:], u)
        x = np.concatenate((np.bitwise_xor(x1, x2), x2))
        out_llr = np.concatenate((out_llr1, out_llr2))
    return u, out_llr, x

def polar_decode_sc(in_llr, F, u):
    N = len(in_llr)
    if N == 2:
        L_1 = f(in_llr[0], in_llr[1])
        u_1 = check(L_1, F[0])
        L_2 = g(in_llr[0], in_llr[1], u_1)
        u_2 = check(L_2, F[1])
        u.append(u_1)
        u.append(u_2)
        return in_llr, u, np.array([u_1^u_2, u_2])
    else:
        
        out_llr1, u1, x1 = polar_decode_sc(f_N(in_llr), F[:int(N/2)], u)
        out_llr2, u, x2 = polar_decode_sc(g_N(in_llr, x1), F[int(N/2):], u1)
        x = np.concatenate((x1^x2, x2))
        out_llr = np.concatenate((out_llr1, out_llr2))
    return out_llr, u, x

def polar_decode_sclist(llrs, L, F):
    list_decode = StructureForList(L, llrs)
    list_decode = polar_decode_sclist_recursive(F, list_decode)
    u_sclist = list_decode.metrics_sequences[0][1]
    return u_sclist

def polar_decode_sclist_crc(llrs, L, F, poly):
    list_decode = StructureForList(L, llrs)
    list_decode = polar_decode_sclist_recursive(F, list_decode)
    u_sclist = list_decode.metrics_sequences[0][1]
    for i, elem in enumerate(list_decode.metrics_sequences):
            if check_crc(np.array(elem[1])[F==0], poly) == 0:
                u_sclist = elem[1]
                break
    return u_sclist

def polar_decode_sclist_recursive(F, list_decode):
    N = len(F)
    if N == 2:
        L = list_decode.length_list()
        list_decode.update_llr('left')
        for i in range(L):
            list_decode.add_left(F[0])
        list_decode.update()
        list_decode.pop_llr()
        
        L = list_decode.length_list()
        list_decode.update_llr('right', depth=1)
        for i in  range(L):
            list_decode.add_right(F[1])
        list_decode.update()
        list_decode.pop_llr()
        
    else:
        list_decode.update_llr('left')
        list_decode = polar_decode_sclist_recursive(F[:int(N/2)], list_decode)
        list_decode.pop_llr()
        
        list_decode.update_llr('right', depth=N//2)
        list_decode = polar_decode_sclist_recursive(F[int(N/2):], list_decode)
        list_decode.pop_llr()
      
        list_decode.update_x(N//2)
    return list_decode
        
        
def polar_decode_scs(n, L, D, in_llr, F):
    stack_decode = StructureForStack(L, D, in_llr, 2**n, F)
    for i in range(n):
        stack_decode.update_llr_top('left')
    
    stack_decode.add()
    
    path_length = stack_decode.get_pathL()
    
    while path_length < 2**n:        
        stack_decode.calc_llr()
        stack_decode.add()
        path_length = stack_decode.get_pathL()
    
    for i in range(n):
        stack_decode.update_x(2**i)
        
    return stack_decode


def polar_decode_scs_crc(n, L, D, in_llr, F, poly):
    current_layer = n
    stack_decode = StructureForStack(L, D, in_llr, 2**n, F)
    #print(np.array(stack_decode.metrics_sequences)[:, :3], 'array len', stack_decode.length_counter, 'stack len', stack_decode.stack_length, 'max_metric', stack_decode.max_metric)
    for i in range(n):
        stack_decode.update_llr_top('left')
    
    stack_decode.add()
    #print(np.array(stack_decode.metrics_sequences)[:, :3], 'array len', stack_decode.length_counter, 'stack len', stack_decode.stack_length, 'max_metric', stack_decode.max_metric)
    
    path_length = stack_decode.get_pathL()
    
    while path_length < 2**n:        
        stack_decode.calc_llr()
        stack_decode.add()
        #print(np.array(stack_decode.metrics_sequences)[:, :3], 'array len', stack_decode.length_counter, 'stack len', stack_decode.stack_length, 'max_metric', stack_decode.max_metric)
        path_length = stack_decode.get_pathL()

    check = check_crc(np.array(stack_decode.metrics_sequences[0][2])[F == 0], poly)
    #print('check', check)
    best_metric_list = stack_decode.metrics_sequences[0]
    ind = 0
    while ind < L and check != 0:
        stack_decode.metrics_sequences.pop(0)
        stack_decode.stack_length -= 1
        if stack_decode.stack_length == 0:
            break
        path_length = stack_decode.get_pathL()
        while path_length < 2**n:        
            stack_decode.calc_llr()
            stack_decode.add()
            #print('after first', np.array(stack_decode.metrics_sequences)[:, :3], 'array len', stack_decode.length_counter, 'stack len', stack_decode.stack_length, 'max_metric', stack_decode.max_metric)
            path_length = stack_decode.get_pathL()
        check = check_crc(np.array(stack_decode.metrics_sequences[0][2])[F == 0], poly)
        #print('check', check)
        ind += 1
    if check != 0:
        stack_decode.metrics_sequences.insert(0, best_metric_list)
        
    for i in range(n):
        stack_decode.update_x(2**i)
    u_scstack = stack_decode.metrics_sequences[0][2]
    #print('final', np.array(stack_decode.metrics_sequences)[:, :3])
    return u_scstack