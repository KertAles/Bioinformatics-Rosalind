# -*- coding: utf-8 -*-
"""
Created on Mon Dec 27 14:52:42 2021

@author: Kert PC
"""


# Computing GC Content


from Bio import Entrez, SeqIO
import sys; sys.path.append(".")

Entrez.email = "ak8754@student.uni-lj.si"


data_read = SeqIO.parse('./rosalind_gc.fasta', 'fasta')

data = {}

for dat in data_read :
    data[dat.id] = {'seq' : str(dat.seq)}
    
maxprop = 0
maxseq = ''
    
for seque in data :
    seq = data[seque]['seq']
    leng = len(seq)
    
    gc_cnt = 0
    for ch in seq :
        if ch == 'G' or ch == 'C' :
            gc_cnt += 1
    
    prop = gc_cnt / leng
    
    if maxprop < prop :
        maxprop = prop
        maxseq = seque
    
    
print(maxseq)
print(maxprop * 100)
