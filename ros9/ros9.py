# -*- coding: utf-8 -*-
"""
Created on Mon Dec 27 15:07:05 2021

@author: Kert PC
"""

from Bio import Entrez, SeqIO
import sys; sys.path.append(".")
from collections import defaultdict
import numpy as np

Entrez.email = "ak8754@student.uni-lj.si"


f = open("./rosalind_ba10c.txt", "r")
x = f.readline()[:-1]
f.readline()
str_alph = " ".join(f.readline().split()).split(' ')
f.readline()
path_alph = " ".join(f.readline().split()).split(' ')
f.readline()
f.readline()

trans = {}
probs = {}

prob_states = {}
prob_states_tmp = {}

for path_ch in path_alph :
    path_split = " ".join(f.readline().split()).split(' ')
    trans[path_split[0]] = {}
    prob_states[path_split[0]] = 0
    prob_states_tmp[path_split[0]] = 0
    
    for idx, path_val in enumerate(path_split[1:]) :
        trans[path_split[0]][path_alph[idx]] = float(path_val) 

f.readline()
f.readline()

for tab_line in f:
    tab_split = " ".join(tab_line.split()).split(' ')
    probs[tab_split[0]] = {}
    
    for idx, tab_val in enumerate(tab_split[1:]) :
        probs[tab_split[0]][str_alph[idx]] = float(tab_val) 


t_1 = np.zeros((len(path_alph), len(x)))
t_2 = np.zeros((len(path_alph), len(x)))

for idx, state in enumerate(path_alph) :
    t_1[idx, 0] = (1.0 / len(path_alph)) * probs[state][x[0]]
        
for idx_1, emit in enumerate(x[1:]) :
    for idx_2, state in enumerate(path_alph) :
        t_1[idx_2, idx_1+1] = max([t_1[idx_3, idx_1] * probs[state][emit] * trans[stat][state] for idx_3, stat in enumerate(path_alph)])
        t_2[idx_2, idx_1+1] = np.argmax([t_1[idx_3, idx_1] * probs[state][emit] * trans[stat][state] for idx_3, stat in enumerate(path_alph)])

j = len(x) - 1

z_t = np.argmax(t_1[:, -1])

ret_str = path_alph[int(z_t)]

while j > 0 :
    z_t = t_2[int(z_t), j]
    ret_str = path_alph[int(z_t)] + ret_str
    j -= 1


with open('res9.txt', 'w') as f:
    f.write(ret_str)

"""
for state in path_alph :
    prob_states[state] = (1.0 / len(path_alph)) * probs[state][x[0]]

with open('res10.txt', 'w') as f:
    f.write("   ".join(path_alph) + '\n')
    for emit in x[1:] :
        sum_probs = 0
        for state_1 in path_alph :
            max_prob = -1
            
            for state_2 in path_alph :
                nu_prob = prob_states[state_2] * trans[state_2][state_1] 

                if nu_prob > max_prob :
                    max_prob = nu_prob
             
            max_prob *= probs[state_1][emit]
            
            prob_states_tmp[state_1] = max_prob
            sum_probs += max_prob
            
        prob_states = prob_states_tmp

        for state in path_alph :
            f.write(str(prob_states[state] / sum_probs) + '  ')
        f.write('\n')
"""
        
        


    


    
