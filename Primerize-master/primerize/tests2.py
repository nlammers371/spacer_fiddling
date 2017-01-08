# -*- coding: utf-8 -*-
"""
Created on Fri Jan 06 10:21:30 2017

@author: Nicholas
"""


sequence = "ggttacccggtacCCGACTGCGTTACTATGGCGACTAtaactgggacaCACATCAAACAC"

N = len(sequence)
m = 20

subset_str1 = []
subset_str2 = []

for i in xrange(N):
    start_pos = max(i - m, 0) - 1
    if start_pos == -1:
        subset_str1.append(sequence[i::-1])
    else:
        subset_str1.append(sequence[i:start_pos:-1])

# match to reverse complement of sequence
sequence_rc = sequence
for i in xrange(N):
    end_pos = N - i - 1
    start_pos = max(end_pos - m - 1, 0) - 1
    if start_pos == -1:
        subset_str2.append(sequence_rc[end_pos::-1])
    else:
        subset_str2.append(sequence_rc[end_pos:start_pos:-1])