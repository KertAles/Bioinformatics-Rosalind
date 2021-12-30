# -*- coding: utf-8 -*-
"""
Created on Mon Dec 27 15:07:05 2021

@author: Kert PC
"""

# Global Alignment with Constant Gap Penalty

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

def global_alignment(seq1, seq2, scoring_function, open_gap_penalty=-5, extend_gap_penalty=0):
    
    indel_char = '*'
    
    seq1 = indel_char + seq1
    seq2 = indel_char + seq2
    
    upp_table = defaultdict(int)
    mid_table = defaultdict(int)
    low_table = defaultdict(int)
    prev = {}
    
    mid_table[0, 0] = 0
    mid_table[0, 1] = open_gap_penalty
    mid_table[1, 0] = open_gap_penalty
    
    upp_table[0, 1] = open_gap_penalty
    low_table[1, 0] = open_gap_penalty
    
    prev[0, 1] = (0, 0)
    prev[1, 0] = (0, 0)
    
    for i in range(2, len(seq1)) :
        mid_table[i, 0] = mid_table[i-1, 0] + extend_gap_penalty
        low_table[i, 0] = low_table[i-1, 0] + extend_gap_penalty
        upp_table[i, 0] = -1000000000000
        prev[i, 0] = (i-1, 0)
    for j in range(2, len(seq2)):
        mid_table[0, j] = mid_table[0, j-1] + extend_gap_penalty
        upp_table[0, j] = upp_table[0, j-1] + extend_gap_penalty
        low_table[0, j] = -1000000000000
        prev[0, j] = (0, j-1)


    for i, si in enumerate(seq1):
        for j, tj in enumerate(seq2):
            if i > 0 and j > 0 :
                low_table[i, j] = max(
                    low_table[i-1, j] + extend_gap_penalty,
                    mid_table[i-1, j] + open_gap_penalty
                )
                upp_table[i, j] = max(
                    upp_table[i, j-1] + extend_gap_penalty,
                    mid_table[i, j-1] + open_gap_penalty
                )
                
                mid_table[i, j], prev[i, j] = max(
                    (mid_table[i-1, j-1] + scoring_function(si, tj), (i-1, j-1)),
                    (upp_table[i, j], (i, j-1)),
                    (low_table[i, j], (i-1, j))
                    )
        
              
    i, j = len(seq1) - 1, len(seq2) - 1
    
    final_score = mid_table[i,j]

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

    return align1, align2, final_score, mid_table



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
s, t, sc, tab = global_alignment(sequences[0], sequences[1], blosum_fun, open_gap_penalty=gap_penalty)

dist = sc

with open('res4.txt', 'w') as f:
    f.write(str(int(dist)))



    
