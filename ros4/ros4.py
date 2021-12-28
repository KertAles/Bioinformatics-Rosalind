# -*- coding: utf-8 -*-
"""
Created on Mon Dec 27 15:07:05 2021

@author: Kert PC
"""

from Bio import Entrez, SeqIO
import sys; sys.path.append(".")
from collections import defaultdict
from Bio.Align import substitution_matrices


from Bio import pairwise2
from Bio.pairwise2 import format_alignment


Entrez.email = "ak8754@student.uni-lj.si"

def hamming_distance(seq1, seq2) :
    dist = 0
    
    for i in range(len(seq1)) :
        if seq1[i] != seq2[i] :
            dist += 1
            
    return dist

def global_alignment(seq1, seq2, scoring_function, constant_gap = True):
    
    table = defaultdict(int)
    prev = {}
    
    indel_char = '*'
    
    seq1 = indel_char + seq1
    seq2 = indel_char + seq2
    
    for i, si in enumerate(seq1):
        for j, tj in enumerate(seq2):
            if i > 0 and j > 0: 
                leftscore = 0
                upscore = 0
                
                if prev[i-1, j] != (i-2, j) :
                    leftscore = scoring_function(si, indel_char)
                    
                if prev[i, j-1] != (i, j-2) :
                    upscore = scoring_function(indel_char, tj)
                
                table[i, j], prev[i, j] = max(
                (table[i-1, j-1] + scoring_function(si, tj), (i-1, j-1)),
                (table[i-1, j] + leftscore, (i-1, j)),
                (table[i, j-1] + upscore, (i, j-1)),
                )
        
                
            elif i == 0 and j > 0:
                table[i, j], prev[i, j] = scoring_function(indel_char, tj), (i, j-1)

            elif i > 0 and j == 0:
                table[i, j], prev[i, j] = scoring_function(si, indel_char), (i-1, j)
    
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

    return align1, align2, final_score, table



data_read = SeqIO.parse('./rosalind_gcon.fasta', 'fasta')
sequences = []
for dat in data_read :
    sequences.append(str(dat.seq))
    

blosum_mat = substitution_matrices.load("blosum62")
gap_penalty = -5

blosum_mat['*', :] = blosum_mat['*', :] - 1
blosum_mat[: , '*'] = blosum_mat[: , '*'] - 1
blosum_mat['*' , '*'] = 1


def align_mat_fun(matrix):
  return lambda s1, s2 : matrix[s1, s2]

blosum_fun = align_mat_fun(blosum_mat)

print('mine')
s, t, sc, tab = global_alignment(sequences[0], sequences[1], blosum_fun)


def gap_function(x, y): 
    if y == 1:  
        return -5
    else : 
        return 0
print('pairwise2')

alignment = pairwise2.align.globalcc(sequences[0], sequences[1], blosum_fun, pairwise2.affine_penalty(-5, 0), pairwise2.affine_penalty(-5, 0), gap_char='*', one_alignment_only=True, score_only=True)

for alg in alignment :
    print(format_alignment(*alg))

#dist = hamming_distance(s, t)
dist = sc
import numpy as np

tabl = np.zeros((len(sequences[0])+1, len(sequences[1])+1))

for idx in tab :
    tabl[idx[0], idx[1]] = tab[idx]


with open('res4.txt', 'w') as f:
    f.write(str(int(dist)))



    
