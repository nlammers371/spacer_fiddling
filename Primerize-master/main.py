# -*- coding: utf-8 -*-
"""
Created on Mon Jan 09 21:45:42 2017

@author: Nicholas
"""

import os
import datetime
import numpy
import primerize

#DEFINE MODEL PARAMETERS

seq1 = 'ggttacccggtactgcataacaatggaacccgaaccgtaactgggacagatcgaaaagctggcctggtttctcgctgtgtgtgccgtgttaatccgtttgccatcagcgagattattagtcaattgcagttgcagcgtttcgctttcgtcctcgtttcactttcgagttagactttattgcagcatcttgaacaatcgtcgcagtttggtaacacgctgtgccatactttcatttagacggaatcgagggaccctggactataatcgcacaacgagaccgggttgcgaagtcagggcattccgccgatctagccatcgccatcttctgcgggcgtttgtttgtttgtttgctgggattagccaagggcttgacttggaatccaatcccgatccctagcccgatcccaatcccaatcccaatcccttgtccttttcattagaaagtcataaaaacacataataatgatgtcgaagggattagggg'
#spacer set 1 used by DePace group with end spacers removed and replaced by wt sequence
seq2 = 'ggttacccggtacCCGACTGCGTTACTATGGCGACTAtaactgggacaCACATCAAACACCATctggtttCACGGTCGATACTAAGTgttaatccgttAGATGATGgcgagattattagtcaattgcagttgcTTGAGATCTTTTAATTCTGTTGAATGGGATCCTAAACgactttattgcagcatcttgaacaatcgtcgcagtttggtaacacgctCATTACTTAGACGCACCagacggaatcgagggaccctggactataatcgcTCCGAATTaccgggttgcGGTACGTTTCTCATTCTGTGACCGGAGGCCTATGTCTGGTATCATGTAGGAgtttgtttgtttgctgggattagccaagggcttgaAGATGGGtccaatcccgatccctagcccgatcccaatcccaatcccaatccctTATTTACTTAGACTgaaagtcataaaaacacataataTCTTTGAcgaagggattagggg'
#seq1 = 'ggttacccggtactgcataacaatggaacccgaaccgtaactgggacagatcgaaaagctggcctggttt'
#seq2 = 'ggttacccggtacCCGACTGCGTTACTATGGCGACTAtaactgggacaCACATCAAACACCATctggttt'
        
min_tm = 60

#test run or results run? 
run = 'tests'
#run = 'results'
#RUN PRIMER DESIGN CODE
job_1d = primerize.Primerize_1D.design(seq1, seq2, MIN_TM=min_tm, NUM_PRIMERS=None, MIN_LENGTH=15, MAX_LENGTH=60, prefix='P4P6_2HP')


if (job_1d.is_success == True):
    #OUTPUT PARAMETERS
    if (run=='results'):
        folderpath = 'Results\\orders'
    elif (run == 'tests'):     
        folderpath = 'Results\\tests'
        
    mydir = os.path.join(os.getcwd(),folderpath, datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
   
    
    if not os.path.exists(mydir):
        os.makedirs(mydir)
        
    #store Inputs
    Inputs = open(os.path.join(mydir,'inputs.txt'), 'w')
    
    InputList = ['Seq1: '+ seq1,'Seq2: '+ seq2, 'Min_Tm: ' + str(min_tm),  job_1d._params]
    
    for item in InputList:
      Inputs.write("%s\n" % item)
    
    Inputs.close()
    
    #save sequences
    Sequences = open(os.path.join(mydir,'sequences.txt'), 'w')
    
    Sequences.write("%s" % job_1d.sequence)
    
    Sequences.close()
    
    
    #save list of primer sequences
    Primers = open(os.path.join(mydir,'primers.txt'), 'w')
    
    Primers.write("%s" % job_1d.primer_set)
    
    Primers.close()
    
    #save corresponding permutation key
    Primer_Permutations = open(os.path.join(mydir,'primer_permutations.txt'), 'w')
    
    Primer_Permutations.write("%s" % job_1d.permutations)
    
    Primer_Permutations.close()
    
    #save corresponding swap key
    Swap_Segments = open(os.path.join(mydir,'swap_segments.txt'), 'w')
    
    Swap_Segments.write("%s" % job_1d.swap_segments)
    
    Swap_Segments.close()
    
    #save array of segment indices
    Segment_Index = open(os.path.join(mydir,'segment_index.txt'), 'w')
    
    Segment_Index.write("%s" % numpy.ndarray.tolist(job_1d.index_array))
    
    Segment_Index.close()
    
    assembly = job_1d._data['assembly']
    
    
    #save primer indices
    Primer_Indices = open(os.path.join(mydir,'primer_indices.txt'), 'w')
    
    Primer_Indices.write("%s" % numpy.ndarray.tolist(assembly.primers))
    
    Primer_Indices.close()