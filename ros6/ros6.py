# -*- coding: utf-8 -*-
"""
Created on Mon Dec 27 15:07:05 2021

@author: Kert PC
"""

# Finding the Longest Multiple Repeat

from Bio import Entrez, SeqIO
import sys; sys.path.append(".")
from collections import defaultdict

Entrez.email = "ak8754@student.uni-lj.si"

def find_largest_substring(graph, parent, k, string) :
    all_children = 0
    all_substrings = []
    
    if parent in graph :
        for child in graph[parent] :
            substrs, num_of_children = find_largest_substring(graph, child['node'], k, string)
            
            if num_of_children >= k :
                for substr in substrs :
                    all_substrings.append(string[child['pos']:child['pos']+child['len']] + substr)
                
                all_substrings.append(string[child['pos']:child['pos']+child['len']])
                    
            all_children += num_of_children
    else :
        all_children += 1

    #print(all_children)
    if all_children < k :
        all_substrings = []
        
    #print(all_substrings)
        
    return (all_substrings, all_children)

f = open("./rosalind_lrep.txt", "r")
string = f.readline()
k = int(f.readline())

suff_graph = defaultdict(list)

for node_line in f:
    node_split = node_line.split(' ')
    suff_graph[node_split[0]].append({'node' : node_split[1], 'pos' : int(node_split[2]) - 1, 'len' : int(node_split[3])})

substrings, bah = find_largest_substring(suff_graph, 'node1', k, string)

maxlen = -1
max_sub = ''
for substr in substrings :
    c_len = len(substr)
    
    if c_len > maxlen :
        maxlen = c_len
        max_sub = substr


with open('res6.txt', 'w') as f:
    f.write(max_sub)


    
