# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 21:27:52 2017

@author: Nicholas
"""

import os
import csv
import re
import numpy
import primerize
import math

#datetime for desired project folder
dt = "2017-01-16_12-08-53"
misprime_warn = 3
#User-generated vector indicating for each swap segment whether to use spacer or wildtype
seg_vec = [0,0,0,0,0,0,0,0,0,0,0,0,1]

#test run or results run?
#run = "tests"
run = "results"

if (run=='results'):
    folderpath = 'Results\\orders'
elif (run == 'tests'):     
    folderpath = 'Results\\tests'
    
inpath = os.path.join(os.getcwd(),folderpath,dt)

assembly_string = ''.join(str(seg_vec))
assembly_string = re.sub('\s','',assembly_string)
out_string = 'assembly_' + assembly_string
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

#generate sequence with segments specified by seg_vec
sequence = ''
swap_ct = 0
for i in xrange(len(index_array[1,:])):
    if (index_array[2,i] < 0):
        sequence += seq1[index_array[0,i]:index_array[1,i]]
    else:
        if (seg_vec[swap_ct]==0):
            sequence += seq1[index_array[0,i]:index_array[1,i]]
        elif (seg_vec[swap_ct]==1):
            sequence += seq2[index_array[0,i]:index_array[1,i]]
        swap_ct += 1

#filter for only vectors that will be used in this assembly
assembly_primers = []
for i in xrange(len(primers)):
    primer_list = primers[i]
    perm_strings = perm_key[i]
    swap_segments = swap_key[i]
    #if primer contains no swap segments include
    if (swap_segments[0] == -1):
          assembly_primers.append(primer_list[0])
    else:
        match_str = ''
        for j in xrange(len(swap_segments)):
            match_str += str(seg_vec[swap_segments[j]])  
            
        for k in xrange(len(primer_list)):
            perm_str = perm_strings[k]
                
            if (perm_str == match_str):
                assembly_primers.append(primer_list[k])
        


(match_list, misprime_score_list, match_id_list) = primerize.misprime._check_misprime_post(sequence, assembly_primers, misprime_warn, primer_key)

#output assembly list
rows = zip(assembly_primers, match_list, misprime_score_list, match_id_list)   
with open(os.path.join(inpath,out_string + '.csv'), 'wb') as f:
    writer = csv.writer(f, delimiter = ',')
    for val in rows:
        writer.writerow(val)


fwd = '-'*len(sequence)
rev = '-'*len(sequence)

for i in xrange(len(primer_key[1,:])):
    primer = assembly_primers[i]
    if (primer_key[2,i]==1):
        fwd = fwd[0:primer_key[0,i]] + primer + fwd[primer_key[1,i]:]
    elif (primer_key[2,i]==-1):
        rev = rev[0:primer_key[0,i]] + primer[::-1] + rev[primer_key[1,i]:]
        
#print assembly
p_len = 150;
p_iter = int(math.ceil(len(sequence)/p_len))

assembly = open(os.path.join(inpath,out_string + '.txt'), 'w')

for i in xrange(p_iter):
    #save list of primer sequences   
    assembly.write("%s \n" % fwd[(i*p_len):(min((i+1)*p_len,len(sequence)))])
    assembly.write("%s \n" % rev[(i*p_len):(min((i+1)*p_len,len(sequence)))])
    assembly.write("\n")    
    
assembly.close()
