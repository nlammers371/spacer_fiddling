# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy
import itertools

def complement(sequence):
    """Convert a DNA sequence input to its complement strand.

    Args:
        sequence: ``str``: Input DNA sequence.

    Returns:
        ``str``: String of complement DNA strand.

    Raises:
        ValueError: For illegal **sequence**.
    """

    rc_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'U': 'A'}
    try:
        sequence = map(lambda x: rc_dict[x], list(sequence))
    except KeyError:
        raise ValueError('\033[41mERROR\033[0m: Illegal sequence value \033[95m%s\033[0m for \033[94mcomplement()\033[0m.\n' % sequence)

    return ''.join(sequence)


def reverse(sequence):
    """Convert a DNA sequence input to its reverse order.

    Args:
        sequence: ``str``: Input DNA sequence.

    returns:
        ``str``: String of reverse DNA strand.
    """

    return sequence[::-1]


def reverse_complement(sequence):
    """Convert a DNA sequence input to its reverse complement strand.

    Args:
        sequence: ``str``: Input DNA sequence.

    Returns:
        ``str``: String of reverse complement DNA strand.
    """

    return complement(reverse(sequence))


""   
"NL Addition: enumerate all possible seuences given 2 sequences"
""

#sequence1 = "ggttacccggtacCCGACTGCGTTACTATGGCGACTAtaactgggacaCACATCAAACAC"
sequence1 = "ggttacccggtactgcataacaatggaacccgaaccgtaactgggacagatcgaaaagct"
sequence2 = "ggttacccggtactgcataacaatggaacccgaaccgtaactgggacagatcgaaaagct"
sequence1 = sequence1.upper()
sequence2 = sequence2.upper()

m=6
l=20
 
subset_str = [] #list of all possible strings
bpID = [] #track starting position of each string
segmentID = [] #track segment number and segment version (list of lists)
strandID = [] #F=0, R=1

for strand in xrange(2):
    if (strand == 0):
        seq1 = sequence1
        seq2 = sequence2
    elif (strand == 1):
        seq1 = reverse_complement(sequence1)
        seq2 = reverse_complement(sequence2)
    
    N = len(seq1)
       
    #ref lists
    cons_list = [] #1 if mismatched region, 0 if matched
    index_list = []
    
    match = 0 #tracks number of consecutive mathing bp
    pmatch = m #previous value of m
    
    for i in xrange(N):
        #if equal, increment by 1, else set to 0
        match = match*(seq1[i]==seq2[i]) + (seq1[i]==seq2[i]) 
        if (match == m):
            #if min conecutive matches reached, backtrack to define start of segement
            index_list.append(i - m + 1)
            cons_list.append(0)
        elif ((match == 0) and (pmatch>=m)):
            #if first mismatch following streak of matches, define new segment
            index_list.append(i) 
            cons_list.append(1)
        pmatch = match
        
    #store start, end, and swap ID of each segment
    index_array = numpy.zeros((3,len(cons_list)), dtype='int')
     
    index_list.append(N)
    
    swap_ID = 1
    keep_ID = -1
    #lists used to store segments forsequence compilation
    swap_list = [] 
    keep_list = []
    
    for i in xrange(len(cons_list)):
        s = index_list[i] 
        e = index_list[i+1]
        index_array[0,i] = s 
        index_array[1,i] = e
        #if mismatched ("swap")region
        if (cons_list[i] == 1):
            index_array[2,i] = swap_ID
            swap_ID += 1
            swap_list.append(seq1[s:e])
            swap_list.append(seq2[s:e])
        
        else:
            index_array[2,i] = keep_ID
            keep_ID += -1 
            keep_list.append(seq1[s:e])
    
    #Iterate through sequence to generate list of all possible strings
    
    output=[]
    
    for j in xrange(1,N+1):
        #j-1 denotes leading edge of string, i trailing
        i = max((j-l),0)
        
        if(strand == 0):
            bp = j-1
        elif(strand == 1):
            bp = N-j

        starts = []
        stops = []
        
        #keep track of position of segments within current sequence
        swap_trk = [] 
        keep_trk = []
        #lookup ids or string list
        swap = [] 
        keep = []
        
        ct = 0
        #find segments that overlap with current sequence
        for k in xrange(len(index_array[1,])):
            
            if ((index_array[0,k] < j) and (index_array[1,k] > i)):
                
                starts.append(index_array[0,k])
                stops.append(index_array[1,k])
    
                if (index_array[2,k] > 0):
                    swap_trk.append(ct) 
                    swap.append((index_array[2,k]-1))
                else:
                    keep_trk.append(ct)
                    keep.append((-index_array[2,k]-1))
                ct += 1
        #if no swap segments, life is simple. There must be only 1 segment involved    
        if (len(swap)==0):
            collapse = keep_list[keep[0]]
            segment = [] #track segment IDs associated with fragment
            segment.append(-9999) #this could be any number. Filler
            
            s_clip = i - min(starts) 
            e_clip = len(collapse) + j - max(stops)
            out = collapse[s_clip:e_clip]
            
            subset_str.append(out[::-1])
            segmentID.append(segment)  
            bpID.append(bp)
            strandID.append(strand)                 
        #otherwise we must enumrate all possible combos         
        else:
            strings = ['']*len(stops)
           
            #fill in segments that are consistent for all segment combinations
            for m in xrange(len(keep)):
                strings[keep_trk[m]] = keep_list[keep[m]]

            #list of enumerated swap combinations in the form of binary strings
            s_rep = ["".join(seq) for seq in itertools.product("01", repeat=len(swap))]
            
            #record each possibility as a separate string
            for m in xrange(len(s_rep)):
                segment = [] #track segment IDs associated with fragment
                binary = s_rep[m] 
                
                for n in xrange(len(swap)):
                    segment.append((swap[n]+1)*(int(binary[n])-(int(binary[n])==0))) # <0 if cons, >0 if swap                    
                    ind = 2*(swap[n]) + int(binary[n])
                    swap_str = swap_list[ind]
                    strings[swap_trk[n]] = swap_str #fill in swap strings
                    
                #concatenate string list    
                collapse = "".join(strings)
                
                s_clip = i - min(starts)
                e_clip = len(collapse) + j - max(stops)
                out = collapse[s_clip:e_clip]

                subset_str.append(out[::-1])
                segmentID.append(segment)                
                bpID.append(bp)
                strandID.append(strand)
    
################################################################
####################Score Nearest Matches#######################

sort_idx = numpy.argsort(subset_str) #track original positions of sorted strings
subset_str.sort() #sort
    
#pre-allocate lists
num_match_mat = numpy.zeros((1, 2*N))-1
misprime_score_mat = numpy.zeros((1, 2*N))-1
best_match_mat = numpy.zeros((1, 2*N))-1


for i in xrange(len(subset_str)):
    seq = subset_str[i] 
    ID = sort_idx[i] #recover original position
    
    segments = segmentID[ID] #segments covered by current seq    
    bp = bpID[ID] + N*strandID[ID] 
 
    #look back
    b = i-1
    match_rev = -1
    misprime_rev = -1
    while (b>=0):
        #check if current comp sequence contains incompatable segments   
        #incompatible segments have IDs of same magnitude but oppoisite sign
        comp_seq = segmentID[sort_idx[b]]        
        comp_seq = [x * -1 for x in comp_seq]        
        if set(segments).isdisjoint(set(comp_seq)):
            count = -1
            misprime_score = 0
            str_1 = seq
            str_2 = subset_str[b]
            #run through strings to see how many consecutive bp's match
            while (count < len(str_1) - 1 and count < len(str_2) - 1): 
                if str_1[count + 1] != str_2[count + 1]:               
                    break
                count += 1
    
                if str_1[count] == 'G' or str_1[count] == 'C':
                    misprime_score += 1.25
                else:
                    misprime_score += 1.0
    
            match_rev = count
            misprime_rev = misprime_score
            break
        else:
            b += -1

    #look forward
    f = i+1
    match_fwd = -1
    misprime_fwd = -1
    while (f<=len(subset_str)-1):
        #check if current comp sequence contains incompatable segments        
        comp_seq = segmentID[sort_idx[f]]        
        comp_seq = [x * -1 for x in comp_seq]        
        if set(segments).isdisjoint(set(comp_seq)):
            count = -1
            misprime_score = 0
            str_1 = seq
            str_2 = subset_str[f]
    
            while (count < len(str_1) - 1 and count < len(str_2) - 1): 
                if str_1[count + 1] != str_2[count + 1]:               
                    break
                count += 1
    
                if str_1[count] == 'G' or str_1[count] == 'C':
                    misprime_score += 1.25
                else:
                    misprime_score += 1.0
    
            match_fwd = count
            misprime_fwd = misprime_score
            break
        else:
            f += 1
            
    if (match_fwd > match_rev):
        match = bpID[sort_idx[f]]
        match_max = match_fwd
        misprime_score_max = misprime_fwd
    else:
        match = bpID[sort_idx[b]]
        match_max = match_rev
        misprime_score_max = misprime_rev
    
    #compare count for current seq to count for position.
    #if higher, take current score
    if (match_max > num_match_mat[0,bp]):
        num_match_mat[0,bp] = match_max
        misprime_score_mat[0,bp] = misprime_score_max
        best_match_mat[0,bp] = match

#split arrays into forward and reverse to match original version output    
num_match_forward = num_match_mat[:,0:N]
misprime_score_forward = misprime_score_mat[:,0:N]
best_match_forward = best_match_mat[:,0:N]

num_match_reverse = num_match_mat[:,N:]
misprime_score_reverse = misprime_score_mat[:,N:]
best_match_reverse = best_match_mat[:,N:]

