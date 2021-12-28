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


f = open("./rosalind_ba10j.txt", "r")
x = f.readline()[:-1]
f.readline()
str_alph = " ".join(f.readline().split()).split(' ')
f.readline()
path_alph = " ".join(f.readline().split()).split(' ')
f.readline()
f.readline()

trans = {}
probs = {}


for path_ch in path_alph :
    path_split = " ".join(f.readline().split()).split(' ')
    trans[path_split[0]] = {}

    for idx, path_val in enumerate(path_split[1:]) :
        trans[path_split[0]][path_alph[idx]] = float(path_val) 

f.readline()
f.readline()

for tab_line in f:
    tab_split = " ".join(tab_line.split()).split(' ')
    probs[tab_split[0]] = {}
    
    for idx, tab_val in enumerate(tab_split[1:]) :
        probs[tab_split[0]][str_alph[idx]] = float(tab_val) 

"""
def forward(seq, state, trans, probs, st_alph):
    summ = 0
    
    if len(seq) > 0 :
        emit = seq[-1]
        seq = seq[:-1]
        
        if len(seq) > 0 :
            summ = 0
            for sta in st_alph :
                summ += forward(seq, sta, trans, probs, st_alph) * trans[sta][state]
        else :
            summ = 1
        
        summ *= probs[state][emit]
    
    return summ
"""

def forward(seq, trans, probs, st_alph):
    emit = seq[0]
    seq = seq[1:]
    
    prob_arr = np.zeros((len(seq)+1, len(st_alph))).astype(float)
    
    for idx, state in enumerate(st_alph) :
        prob_arr[0, idx] = probs[state][emit] * (1.0 / len(st_alph))

    i = 1
    while len(seq) > 0 :
        emit = seq[0]
        
        for idx, state in enumerate(st_alph) :
            for idx2, state_prev in enumerate(st_alph) :
                 prob_arr[i, idx] += prob_arr[i-1, idx2] * trans[state_prev][state]
                 
            prob_arr[i, idx] *= probs[state][emit]
            
        seq = seq[1:]
        i += 1
    
    return prob_arr


def backward(seq, trans, probs, st_alph):
    #emit = seq[-1]
    #seq = seq[:1]
    
    prob_arr = np.zeros((len(seq), len(st_alph))).astype(float)
    
    for idx, state in enumerate(st_alph) :
        prob_arr[-1, idx] = 1

    i = len(seq) - 2
    while len(seq) > 1 :
        emit = seq[-1]
        #summ = 0
        for idx, state in enumerate(st_alph) :
            for idx2, state_next in enumerate(st_alph) :
                prob_arr[i, idx] += trans[state][state_next] * probs[state_next][emit] * prob_arr[i+1, idx2]
            #summ += prob_arr[i, idx]
            
        #prob_arr[i, :] /= summ
        seq = seq[:-1]
        i -= 1
    
    emit = seq[0]
    prob_test = 0
    for idx, state in enumerate(st_alph) :
        prob_test += (1 / len(st_alph)) * probs[state][emit] * prob_arr[0, idx]
    
    print(prob_test)
    
    return prob_arr

"""
def backward(seq, trans, probs, st_alph):
    #emit = seq[-1]
    #seq = seq[:1]
    
    prob_arr = np.zeros((len(seq), len(st_alph))).astype(float)

    for i, ch_i_plus in enumerate(seq[1:][::-1]):
        
        for idx, state in enumerate(st_alph) :
            if i == 0:
                prob_arr[i, idx] = 1 / len(st_alph)
            else :
                prob_arr[i, idx] = np.sum(trans[state][l] * probs[l][ch_i_plus] * prob_arr[i-1, idx] for l in st_alph)
    
    return prob_arr
"""

prob_arr = forward(x, trans, probs, path_alph)
prob_arr_bac = backward(x, trans, probs, path_alph)

prob_x = np.sum(prob_arr[-1, :])
#prob_x_bac = np.sum(prob_arr_bac[0, :])


with open('res10.txt', 'w') as f:
    f.write("\t".join(path_alph) + '\n')
    for idx, ch in enumerate(x) :
        for idx2, sta in enumerate(path_alph) :
            prob_state = (prob_arr[idx][idx2] * prob_arr_bac[idx][idx2]) / prob_x
            
            
            f.write(("%.4f" % round(prob_state, 4)).rstrip('0'))
            if idx2 < len(path_alph)-1 :
                f.write('\t')
        f.write('\n')
        
f.close()
print('Done')
        

    
