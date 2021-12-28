# -*- coding: utf-8 -*-
"""
Created on Mon Dec 27 15:07:05 2021

@author: Kert PC
"""

# Compute the Probability of an Outcome Given a Hidden Path

from Bio import Entrez, SeqIO
import sys; sys.path.append(".")
from collections import defaultdict

Entrez.email = "ak8754@student.uni-lj.si"


f = open("./rosalind_ba10b.txt", "r")
x = f.readline()[:-1]
f.readline()
str_alph = " ".join(f.readline().split()).split(' ')
f.readline()
path = f.readline()[:-1]
f.readline()
path_alph = " ".join(f.readline().split()).split(' ')
f.readline()
f.readline()

probs = {}

for tab_line in f:
    tab_split = " ".join(tab_line.split()).split(' ')
    probs[tab_split[0]] = {}
    
    for idx, tab_val in enumerate(tab_split[1:]) :
        probs[tab_split[0]][str_alph[idx]] = float(tab_val) 

prob = 1.0

for idx, ch in enumerate(x) :
    prob = prob * probs[path[idx]][ch]

result = prob

with open('res7.txt', 'w') as f:
    f.write(str(result))


    
