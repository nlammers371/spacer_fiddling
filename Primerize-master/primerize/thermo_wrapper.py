# -*- coding: utf-8 -*-
"""
Created on Tue Jan 03 13:19:00 2017

@author: Nicholas
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jan 03 13:19:00 2017

@author: Nicholas
"""

from . import thermo
import numpy
import itertools

#seq1 = 'ggttacccggtactgcataacaatggaacccgaaccgtaactgggacagatcgaaaagctggcctggtttctcgctgtgtgtgccgtgttaatccgtttgccatcagcgagattattagtcaattgcagttgcagcgtttcgctttcgtcctcgtttcactttcgagttagactttattgcagcatcttgaacaatcgtcgcagtttggtaacacgctgtgccatactttcatttagacggaatcgagggaccctggactataatcgcacaacgagaccgggttgcgaagtcagggcattccgccgatctagccatcgccatcttctgcgggcgtttgtttgtttgtttgctgggattagccaagggcttgacttggaatccaatcccgatccctagcccgatcccaatcccaatcccaatcccttgtccttttcattagaaagtcataaaaacacataataatgatgtcgaagggattagggg'
#seq2 = 'ggttacccggtacCCGACTGCGTTACTATGGCGACTAtaactgggacaCACATCAAACACCATctggtttCACGGTCGATACTAAGTgttaatccgttAGATGATGgcgagattattagtcaattgcagttgcTTGAGATCTTTTAATTCTGTTGAATGGGATCCTAAACgactttattgcagcatcttgaacaatcgtcgcagtttggtaacacgctCATTACTTAGACGCACCagacggaatcgagggaccctggactataatcgcTCCGAATTaccgggttgcGGTACGTTTCTCATTCTGTGACCGGAGGCCTATGTCTGGTATCATGTAGGAgtttgtttgtttgctgggattagccaagggcttgaAGATGGGtccaatcccgatccctagcccgatcccaatcccaatcccaatccctTATTTACTTAGACTgaaagtcataaaaacacataataTCTTTGAcgaagggattagggg'
#seq1 = 'ggttacccggtactgcataacaatggaac'
#seq2 = 'ggttacccggtacCCGACTGCGTTACTAT'
def _thermo_multi_seq(seq1, seq2, min_length=15, max_length=60, m=6):
    seq1 = seq1.upper()
    seq2 = seq2.upper()
   
    N = len(seq1)
   
   #ref lists
    cons_list = []
    ID_list = []
    index_list = []
    
    ID = 0
    match = 0
    pmatch = m

    for i in xrange(N):
        #if equal, increment by 1, else set to 0
        match = match*(seq1[i]==seq2[i]) + (seq1[i]==seq2[i]) 
        if (match == m):
            #if min conecutive matches reached, backtrack to define start of seq
            index_list.append(i - m + 1)
            ID_list.append(ID)
            cons_list.append(0)
            ID += 1
        elif ((match == 0) and (pmatch>=m)):
            index_list.append(i) 
            ID_list.append(ID)
            cons_list.append(1)
            ID += 1
        pmatch = match
        
    #store start, end, and swap ID of each segment
    index_array = numpy.zeros((3,len(cons_list)), dtype='int')
     
    index_list.append((N))
    
    swap_ID = 1
    keep_ID = -1
    swap_list = []
    keep_list = []
    
    for i in xrange(len(cons_list)):
        s = index_list[i] 
        e = index_list[i+1]
        index_array[0,i] = s 
        index_array[1,i] = e
        if (cons_list[i] == 1):
            index_array[2,i] = swap_ID
            swap_ID += 1
            swap_list.append(seq1[s:e])
            swap_list.append(seq2[s:e])
        else:
            index_array[2,i] = keep_ID
            keep_ID += -1 
            keep_list.append(seq1[s:e])
    
    #create lookupt table for TM values
    #when there are multiple distinct TM scores for a position, take minimum
    tm_precalc = numpy.zeros((N,N),dtype='float') - 99999
    for i in xrange((N-min_length)+1):
        
         for j in xrange((i+min_length),min(N,(i+max_length+1))):
            
            IDs = []
            starts = []
            stops = []
            
            swap_trk = []
            swap = []
            keep_trk = []
            keep = []
            
            ct = 0
            for k in xrange(len(index_array[1,])):
                
                if ((index_array[0,k] < j) and (index_array[1,k] > i)):
                    IDs.append(k)
                    starts.append(index_array[0,k])
                    stops.append(index_array[1,k])
                    if (index_array[2,k] > 0):
                        swap_trk.append(ct) 
                        swap.append((index_array[2,k]-1))
                    else:
                        keep_trk.append(ct)
                        keep.append((-index_array[2,k]-1))
                    ct += 1
                    
            tm_scores = []
            if (len(swap)==0):
                collapse = keep_list[keep[0]]
                
                s_clip = i - min(starts)
                e_clip = len(collapse) + j - max(stops)
                
                out = collapse[s_clip:e_clip]
                tm = thermo.calc_Tm(out)
                tm_scores.append(tm)
            else:
                strings = ['']*len(IDs)
                
                for m in xrange(len(keep)):
                    strings[keep_trk[m]] = keep_list[keep[m]]
    
                s_rep = ["".join(seq) for seq in itertools.product("01", repeat=len(swap))]
                for m in xrange(len(s_rep)):
                    bin = s_rep[m]
                    for n in xrange(len(swap)):
                        ind = 2*(swap[n]) + int(bin[n])
                        swap_str = swap_list[ind]
                        strings[swap_trk[n]] = swap_str
                        
                    collapse = "".join(strings)
                    
                    s_clip = i - min(starts)
                    e_clip = len(collapse) + j - max(stops)
                    
                    out = collapse[s_clip:e_clip]
                    tm = thermo.calc_Tm(out)
                    tm_scores.append(tm)

            tm_precalc[i,(j-1)] = min(tm_scores)
            
    return(tm_precalc)
                
            
