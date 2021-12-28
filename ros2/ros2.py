# -*- coding: utf-8 -*-
"""
Created on Mon Dec 27 15:07:05 2021

@author: Kert PC
"""

from Bio import Entrez, SeqIO
import sys; sys.path.append(".")
from collections import defaultdict

Entrez.email = "ak8754@student.uni-lj.si"

def hamming_distance(seq1, seq2) :
    dist = 0
    
    for i in range(len(seq1)) :
        if seq1[i] != seq2[i] :
            dist += 1
            
    return dist

def global_alignment(seq1, seq2, scoring_function):
    
    table = defaultdict(int)
    prev = {}
    
    indel_char = '-'
    
    seq1 = indel_char + seq1
    seq2 = indel_char + seq2
    
    for i, si in enumerate(seq1):
        for j, tj in enumerate(seq2):
            if i > 0 and j > 0: 
                table[i, j], prev[i, j] = max(
                (table[i-1, j] + scoring_function(si, indel_char), (i-1, j)),
                (table[i, j-1] + scoring_function(indel_char, tj), (i, j-1)),
                (table[i-1, j-1] + scoring_function(si, tj), (i-1, j-1))
                )
            elif i == 0 and j > 0:
                table[i, j], prev[i, j] = table[i, j-1] + scoring_function(indel_char, tj), (i, j-1)

            elif i > 0 and j == 0:
                table[i, j], prev[i, j] = table[i-1, j] + scoring_function(si, indel_char), (i-1, j)
    
    i, j = len(seq1)-1, len(seq2)-1
    
    final_score = table[i,j]

    align1 = ''
    align2 = ''

    while (i, j) != (0,0) :
        if prev[i, j] == (i-1, j-1):
            align1 = seq1[i] + align1
            align2 = seq2[j] + align2
        elif prev[i, j] == (i, j-1):
            align1 = '-' + align1
            align2 = seq2[j] + align2
        elif prev[i, j] == (i-1, j):
            align1 = seq1[i] + align1
            align2 = '-' + align2
            
        i, j = prev[i, j]

    return align1, align2, final_score 



data_read = SeqIO.parse('./rosalind_edta.fasta', 'fasta')

sequences = []

for dat in data_read :
    sequences.append(str(dat.seq))
    
s, t, sc = global_alignment(sequences[0], sequences[1], lambda x, y: [-1, 0][x == y])

#dist = hamming_distance(s, t)
dist = abs(sc)


with open('res2.txt', 'w') as f:
    f.write(str(dist))
    f.write('\n')
    f.write(s)
    f.write('\n')
    f.write(t)
    f.write('\n')



    
