import numpy as np
from utils import check_list, f, f_N, g, g_N


class StructureForList:
    def __init__(self, L, in_llr):
        self.L = L
        self.metrics_sequences = [[0, [], [in_llr], np.array([])]]     
    
    def add_left(self, F):
        metric_old, u_old, llr, x = self.metrics_sequences.pop(0)
        Li = llr[-1][0]
        check_res = check_list(Li, F)
        for bit, metric in check_res:
            u = u_old.copy()
            u.append(bit)
            x1 = np.concatenate((x,  [bit]))
            self.metrics_sequences.append([metric_old+metric, u, llr.copy(), x1])
        
    def add_right(self, F):
        metric_old, u_old, llr, x = self.metrics_sequences.pop(0)
        Li = llr[-1][0]
        check_res = check_list(Li, F)
        for bit, metric in check_res:
            u = u_old.copy()
            bit_prev = u[-1]
            u.append(bit)
            x1 = np.concatenate((x[:-1],  [bit_prev^bit, bit]))
            self.metrics_sequences.append([metric_old+metric, u, llr.copy(), x1])
    
    def update(self):
        self.metrics_sequences = sorted(self.metrics_sequences, key=lambda metric: metric[0])
        self.metrics_sequences = self.metrics_sequences[:self.L]
        
    def length_list(self):
        return len(self.metrics_sequences)
    
    def pop_llr(self):
        for i, elem in enumerate(self.metrics_sequences):
            self.metrics_sequences[i][2].pop()
    
    def update_llr(self, l_r, depth=None):
        length = self.length_list()
        if l_r == 'left':
            for i in range(length):
                llr = self.metrics_sequences[i][2][-1]
                self.metrics_sequences[i][2].append(f_N(llr))
        elif l_r == 'right':
            for i in range(length):
                llr = self.metrics_sequences[i][2][-1]
                self.metrics_sequences[i][2].append(g_N(llr, self.metrics_sequences[i][3][-depth:]))
            
    def update_x(self, depth):
        for i, elem in enumerate(self.metrics_sequences):
            x = elem[3].astype(int)
            self.metrics_sequences[i][3] = np.concatenate((x[:-2*depth],x[-2*depth:-depth]^x[-depth:], x[-depth:]))
            

class StructureForStack:
    def __init__(self, L, D, in_llr, N, F):
        self.L = L
        self.D = D
        self.F = F
        self.metrics_sequences = [[0, 0, [], [in_llr], np.array([])]] 
        self.length_counter = np.zeros(N)
        self.stack_length = 1
        self.max_metric = 0

    
    def add(self):
        metric_old, length_old, u_old, llr, x = self.metrics_sequences.pop(0)
        self.stack_length -= 1
        self.length_counter[length_old] += 1
        if self.length_counter[length_old] == self.L: 
            for i, elem in enumerate(self.metrics_sequences):
                if elem[1] <= length_old:
                    self.metrics_sequences.pop(i)
                    self.stack_length -= 1
        if(self.stack_length != 0):
            self.max_metric = self.metrics_sequences[-1][0]
        else:
            self.max_metric = 0
        Li = llr[-1][0]
        check_res = check_list(Li, self.F[length_old])
        for bit, metric in check_res:
            u = u_old.copy()
            u.append(bit)
            x1 = np.concatenate((x,  [bit]))
            new_metric = metric_old+metric
            if new_metric < self.max_metric:
                for i, elem in enumerate(self.metrics_sequences):
                    if new_metric <= elem[0]:
                        self.metrics_sequences.insert(i, [new_metric, length_old+1, u, llr.copy(), x1])
                        self.stack_length += 1
                        break
            elif self.stack_length <= self.D:
                self.metrics_sequences.append([new_metric, length_old+1, u, llr.copy(), x1])
                self.max_metric = new_metric
                self.stack_length += 1
        while self.stack_length > self.D:
            self.metrics_sequences.pop()
            self.stack_length -= 1
        self.max_metric = self.metrics_sequences[-1][0]
            
    def calc_llr(self):
        pathL = self.get_pathL()
        ind = pathL
        layer = 0
        self.pop_llr_top()
        while ind % 2 == 0:
            self.update_x(2**layer)
            self.pop_llr_top()
            layer += 1
            ind /= 2 
        self.update_llr_top('right', depth = 2**layer)
        for i in range(layer):
            self.update_llr_top('left')
        
        
    def length_stack(self):
        return len(self.metrics_sequences)
    
    def get_pathL(self):
        return self.metrics_sequences[0][1]
    
    def pop_llr_top(self):
        self.metrics_sequences[0][3].pop()
    
    def update_llr_top(self, l_r, depth=None):
        if l_r == 'left':
            llr = self.metrics_sequences[0][3][-1]
            self.metrics_sequences[0][3].append(f_N(llr))
        elif l_r == 'right':
            llr = self.metrics_sequences[0][3][-1]
            self.metrics_sequences[0][3].append(g_N(llr, self.metrics_sequences[0][4][-depth:]))
            
    def update_x(self, depth):
        x = self.metrics_sequences[0][4].copy().astype(int)
        self.metrics_sequences[0][4] = np.concatenate((x[:-2*depth],x[-2*depth:-depth]^x[-depth:], x[-depth:]))

