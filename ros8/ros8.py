# -*- coding: utf-8 -*-
"""
Created on Mon Dec 27 15:07:05 2021

@author: Kert PC
"""

# Compute the Probability of a String Emitted by an HMM

from Bio import Entrez, SeqIO
import sys; sys.path.append(".")
from collections import defaultdict
import numpy as np

Entrez.email = "ak8754@student.uni-lj.si"


f = open("./rosalind_ba10d.txt", "r")
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
    emit = seq[0]
    seq = seq[1:]
    
    if len(seq) > 0 :
        summ = 0
        for sta in st_alph :
            summ += forward(seq, sta, trans, probs, st_alph, prev_probs) * trans[sta][state]
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


prob_arr = forward(x, trans, probs, path_alph)

prob = np.sum(prob_arr[-1, :])

with open('res8.txt', 'w') as f:
    f.write(str(prob))



