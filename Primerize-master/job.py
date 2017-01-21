# -*- coding: utf-8 -*-
"""
Created on Mon Jan 09 22:34:06 2017

@author: Nicholas
"""
import os
import re
import numpy
#import primerize
import csv
#import datetime
#import json
#datetime for desired project folder
dt = "2017-01-16_12-08-53"
   
#For first pass I want to test swap segments  4 and 11 since they are the most highly conserved
seg_vec = [[0],[0],[0],[0,1],[0],[0],[0],[0],[0],[0],[0,1]]


inpath = os.path.join(os.getcwd(),'Results\\orders',dt)

#read in sequences
Sequences = open(os.path.join(inpath,'sequences.txt'), 'r')

seq_str = Sequences.read()
seq_list = eval(seq_str)
seq1 = seq_list[0]
seq2 = seq_list[1]
Sequences.close()

#read in primer list
Primers = open(os.path.join(inpath,'primers.txt'), 'r')

prime_list = Primers.read()
primers = eval(prime_list)
Primers.close()

#read in primer permutation key
Perm_key = open(os.path.join(inpath,'primer_permutations.txt'), 'r')

perm_list = Perm_key.read()
perm_key = eval(perm_list)
Perm_key.close()


#read in primer swap segment key
Swap_key = open(os.path.join(inpath,'swap_segments.txt'), 'r')

swap_list = Swap_key.read()
swap_key = eval(swap_list)
Swap_key.close()

#read in primer index array
Primer_key = open(os.path.join(inpath,'primer_indices.txt'), 'r')

primer_str = Primer_key.read()
primer_key = numpy.asarray(eval(primer_str))
Primer_key.close()

#read in index array
Index_array = open(os.path.join(inpath,'segment_index.txt'), 'r')

index_list = Index_array.read()
index_array = numpy.asarray(eval(index_list))

Index_array.close()


swap_strings = [['']]*len(swap_key)

for i in xrange(len(swap_key)):
    if (max(swap_key[i]) == -1):
        swap_strings[i] = ['0']
        
for i in xrange(len(seg_vec)):
    swap_bin = seg_vec[i]
    for j in xrange(len(swap_key)):
        swap = swap_key[j]
        swap_string = swap_strings[j]
        if (i in swap):
            swap_new = []
            for k in xrange(len(swap_string)):
                
                for l in xrange(len(swap_bin)):
                    swap_new.append(swap_string[k]+str(swap_bin[l]))
                    
            swap_strings[j] = swap_new
            

primer_index = []           
for i in xrange(len(perm_key)):
    permutations = perm_key[i]
    perm_out = []
    for j in xrange(len(permutations)):
       if (permutations[j] in swap_strings[i]):
           perm_out.append(j)
    primer_index.append(perm_out)
    
primer_set = []
primer_id = []
primer_key_start = []
primer_key_stop = []
strand_id = []

for i in xrange(len(primers)):
    primer_list = primers[i]
    ind_list = primer_index[i]
    swap_num_list = swap_key[i]
    swap_str_list = swap_strings[i]
    for j in xrange(len(ind_list)):
        primer_key_start.append(primer_key[0,i])
        primer_key_stop.append(primer_key[1,i]) 
        strand_id.append(primer_key[2,i])
        swap_id = []        
        primer_set.append(primer_list[ind_list[j]])
        swap_str = swap_str_list[j]
        for k in xrange(len(swap_str)):
            swap_id.append((0-int(swap_str[k])+int(swap_str[k]=='0'))*(1+swap_num_list[k]))
        primer_id.append(swap_id)
        
#(match_mat, warn_scores, warn_indices) = primerize.misprime._check_misprime_post(primer_set, primer_id, misprime_warn)

###############################################################################
####Output Results

job_string = ''.join(str(seg_vec))
job_string = re.sub('\s','',job_string)
outpath = os.path.join(inpath, 'Job ' + job_string)

if not os.path.exists(outpath):
    os.makedirs(outpath)
 

rows = zip(primer_set, primer_key_start, primer_key_stop, primer_id, strand_id)   
with open(os.path.join(outpath,'primer_list.csv'), 'wb') as f:
    writer = csv.writer(f, delimiter = ',')
    for val in rows:
        writer.writerow(val)
        
bp_total = len("".join(primer_set))
print(str(len(seq1)) + 'bp in original sequence')
print(str(bp_total) + 'bp in primer set' )