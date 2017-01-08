# -*- coding: utf-8 -*-
"""
Created on Mon Jan 02 21:51:32 2017

@author: Nicholas
"""

import numpy

def segment(seq1,seq2,m=6):
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    
    N = len(seq1)
    
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
    
    #index_list.append(N)
    
    #store start, end, and swap ID of each segment
    index_array = numpy.zeros((3,len(cons_list)), dtype='int')
     
    index_list.append((N))
    
    
    for i in xrange(len(cons_list)):
        
        index_array[0,i] = index_list[i] 
        index_array[1,i] = index_list[i+1]
        index_array[2,i] = cons_list[i+1]        
        
    
    #array to track number of swap segments between any two points
    swap_mat = numpy.zeros((N,N),dtype='float')-9999
    for i in xrange(N-1):
        for j in xrange((i+1),N+1): 
            cons = []
            
            for k in xrange(len(index_array[1,])):          
                if ((index_array[0,k] < j) and (index_array[1,k] > i)):
                    cons.append(index_array[2,k])
                
            swap_mat[i,j-1] = sum(cons)
   
    return(swap_mat)             
        
    
    
    
      