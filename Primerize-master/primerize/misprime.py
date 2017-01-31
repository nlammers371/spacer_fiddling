import numpy
import itertools
from .util import reverse_complement


def _check_misprime(sequence):
    N = len(sequence)
    # length of string subsets
    m = 20
    subset_str = []

    # match to sequence
    for i in xrange(N):
        start_pos = max(i - m, 0) - 1
        if start_pos == -1:
            subset_str.append(sequence[i::-1])
        else:
            subset_str.append(sequence[i:start_pos:-1])

    # match to reverse complement of sequence
    sequence_rc = reverse_complement(sequence)
    for i in xrange(N):
        end_pos = N - i - 1
        start_pos = max(end_pos - m - 1, 0) - 1
        if start_pos == -1:
            subset_str.append(sequence_rc[end_pos::-1])
        else:
            subset_str.append(sequence_rc[end_pos:start_pos:-1])

    sort_idx = numpy.argsort(subset_str)
    subset_str.sort()

    # how close is match to neighbor?

    match_next = numpy.zeros((1, 2 * N - 1))
    misprime_score_next = numpy.zeros((1, 2 * N - 1))
    for i in xrange(2 * N - 1):
        count = -1
        misprime_score = 0
        str_1 = subset_str[i]
        str_2 = subset_str[i + 1]

        while (count < len(str_1) - 1 and count < len(str_2) - 1): 
            if str_1[count + 1] != str_2[count + 1]:               
                break
            count += 1

            if str_1[count] == 'G' or str_1[count] == 'C':
                misprime_score += 2
            else:
                misprime_score += 1.0

        match_next[0, i] = count
        misprime_score_next[0, i] = misprime_score

    match_max = numpy.zeros((1, 2 * N))
    best_match = numpy.zeros((1, 2 * N), dtype=numpy.int16)
    misprime_score_max = numpy.zeros((1, 2 * N))

    # compare both neighbors
    match_max[0, 0] = match_next[0, 0]
    best_match[0, 0] = 1
    misprime_score_max[0, 0] = misprime_score_next[0, 0]

    match_max[0, 2 * N - 1] = match_next[0, 2 * N - 2]
    best_match[0, 2 * N - 1] = 2 * N - 2
    misprime_score_max[0, 2 * N - 1] = misprime_score_next[0, 2 * N - 2]

    for i in xrange(1, 2 * N - 1):
        if (match_next[0, i - 1] > match_next[0, i]):
            best_match[0, i] = i - 1
            match_max[0, i] = match_next[0, i - 1]
            misprime_score_max[0, i] = misprime_score_next[0, i - 1]
        else:
            best_match[0, i] = i + 1
            match_max[0, i] = match_next[0, i]
            misprime_score_max[0, i] = misprime_score_next[0, i]

    num_match_foward = numpy.zeros((1, N))
    misprime_score_forward = numpy.zeros((1, N))
    best_match_forward = numpy.zeros((1, N))
    num_match_reverse = numpy.zeros((1, N))
    misprime_score_reverse = numpy.zeros((1, N))
    best_match_reverse = numpy.zeros((1, N))

    for i in xrange(2 * N):
        if (sort_idx[i] <= N - 1):
            num_match_foward[0, sort_idx[i]] = match_max[0, i]
            misprime_score_forward[0, sort_idx[i]] = misprime_score_max[0, i]
            best_match_forward[0, sort_idx[i]] = (sort_idx[best_match[0, i]] - 1) % N + 1
        else:
            num_match_reverse[0, sort_idx[i] - N] = match_max[0, i]
            misprime_score_reverse[0, sort_idx[i] - N] = misprime_score_max[0, i]
            best_match_reverse[0, sort_idx[i] - N] = (sort_idx[best_match[0, i]] - 1) % N + 1

    return (num_match_foward, num_match_reverse, best_match_forward, best_match_reverse, misprime_score_forward, misprime_score_reverse)

""
"NL Addition: multi-seq version of misprime function"
""


def _check_misprime_multi_archive(sequence1, sequence2, m=6, l=20):
    
    sequence1 = sequence1.upper()
    sequence2 = sequence2.upper()        
        
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
                        swap_str = swap_list[ind] #wt if zero, swap spacer if 1
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
     
        #look back until a compatible comparison seq is found
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
                        misprime_score += 2.0 #NL: Studies indicate that energy associated with GC bonds is about twice that of AT  
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
                        misprime_score += 2.0 #NL: Studies indicate that energy associated with GC bonds is about twice that of AT  
                    else:
                        misprime_score += 1.0
        
                match_fwd = count
                misprime_fwd = misprime_score
                break
            else:
                f += 1
                
        if (match_rev > match_fwd):
            match = bpID[sort_idx[b]]
            match_max = match_rev
            misprime_score_max = misprime_rev
        else:
            match = bpID[sort_idx[f]]
            match_max = match_fwd
            misprime_score_max = misprime_fwd
        
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
    
    return (num_match_forward, num_match_reverse, best_match_forward, best_match_reverse, misprime_score_forward, misprime_score_reverse, subset_str)


def _check_misprime_multi(sequence1, sequence2, m=6, l=20):
    
    sequence1 = sequence1.upper()
    sequence2 = sequence2.upper()        
        
    subset_str = [] #list of all possible sequences
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
        
        N_BP = len(seq1)
           
        #ref lists
        cons_list = [] #1 if mismatched region, 0 if matched
        index_list = [] #starting bp of each new segments
        
        match = 0 #tracks number of consecutive matching bp
        pmatch = m #previous value of m
        #with current formailism, it is impossible to have a "match" segment
        #<m bp in lengths but it is possible to have a shorter mismatched segment
        for i in xrange(N_BP):
            #if equal, increment by 1, else set to 0
            match = match*(seq1[i]==seq2[i]) + (seq1[i]==seq2[i]) 
            if (match == m):
                #if min conecutive matches reached, backtrack to define start of segement
                index_list.append(i - m + 1)
                cons_list.append(0)
            elif ((match == 0) and (pmatch >= m)):
                #if first mismatch following streak of matches, define new segment
                index_list.append(i) 
                cons_list.append(1)
            pmatch = match
            
        #store start, end, and segment ID and swap ID of each segment
        index_array = numpy.zeros((4,len(cons_list)), dtype='int')
         
        index_list.append((N_BP))
        #account for fact that we are working with RC. Segment numbering will need 
        #to be consistent for subsequent misprime scoring        
    
        swap_ID = 1
        #store sequence strings
        swap_list = []        
        
        for i in xrange(len(cons_list)):
            s = index_list[i] 
            e = index_list[i+1]
            index_array[0,i] = s 
            index_array[1,i] = e
            if (cons_list[i] == 1):
                index_array[2,i] = swap_ID #convention: swap segments > 0, keep < 0
                swap_list.append([seq1[s:e],seq2[s:e]])
            else:
                index_array[2,i] = -swap_ID
                swap_list.append(seq1[s:e])
            index_array[3,i] = swap_ID
            swap_ID += 1
            
        #if on reverse strand we need an index var to recover correct segment ID
        segment_ind = strand*(1+len(cons_list))
            
        #Iterate through sequence to generate list of all possible strings
        for j in xrange(1,N_BP+1):
            #j-1 denotes leading edge of string, i trailing
            i = max((j-l),0)
            
            if(strand == 0):
                bp = j-1
            elif(strand == 1):
                bp = N_BP-j
    
           #filter for relevant segments
            index_filt = index_array[:, (index_array[0,:] < j) & (index_array[1,:] > i) ]
            ind_start = min(index_filt[3,:]) 
            
            #break up array for simplicity
            starts = index_filt[0,:]
            stops = index_filt[1,:]
            
            
            #list to store segment strings            
            strings = ['']*len(starts)  
            
            keep_indices = index_filt[3,(index_filt[2,:] < 0)]-ind_start
            swap_indices = index_filt[3,(index_filt[2,:] > 0)]-ind_start
 
            for k in keep_indices:
                strings[k] = swap_list[k+ind_start-1]
                     
            #if no swap segments, life is simple. There must be only 1 segment involved                 
            if (len(swap_indices)==0):
                collapse = "".join(strings)
                s_clip = i - min(starts)
                e_clip = len(collapse) + j - max(stops)
                segment = [-9999] #this is filler. could be any number
                                
                out = collapse[s_clip:e_clip]
                
                subset_str.append(out[::-1])
                segmentID.append(segment)  
                bpID.append(bp)
                strandID.append(strand)   
                
            #otherwise we must enumrate all possible combos         
            elif (len(swap_indices)>0):                       
                #for len(swap) swap segments, enumerate all possible unique combinations
                s_rep = ["".join(seq) for seq in itertools.product("01", repeat=len(swap_indices))]
                for n in xrange(len(s_rep)):
                    segment = []
                    bin = s_rep[n]
                    for o in xrange(len(bin)):
                        seg_ID = abs(segment_ind - (swap_indices[o] + ind_start))
                        segment.append(seg_ID*(int(bin[o])-(int(bin[o])==0))) # <0 if cons, >0 if swap
                        swap_strings = swap_list[swap_indices[o] + ind_start - 1]
                        swap_str = swap_strings[int(bin[o])]
                        strings[swap_indices[o]] = swap_str
                
                #concatenate string list    
                    collapse = "".join(strings)
                    
                    s_clip = i - min(starts)
                    e_clip = len(collapse) + j - max(stops)
                    out = collapse[s_clip:e_clip]
    
                    subset_str.append(out[::-1])
                    segmentID.append(segment)                
                    bpID.append(bp)
                    strandID.append(strand)
            else:
                print('No swap or keep segments found for index ' + str(i) + ',' + str(j) + ' on strand ' + str(strand))
    
    ###############################################################
    ####################Score Nearest Matches#######################
    
    sort_idx = numpy.argsort(subset_str) #track original positions of sorted strings
    subset_str.sort() #sort sequences fragments
        
    #pre-allocate lists
    num_match_mat = numpy.zeros((1, 2*N_BP))-1
    misprime_score_mat = numpy.zeros((1, 2*N_BP))-1
    best_match_mat = numpy.zeros((1, 2*N_BP))-1    
    
    for i in xrange(len(subset_str)):
        seq = subset_str[i] 
        ID = sort_idx[i] #recover original position
        
        segments = segmentID[ID] #segments covered by current seq    
        bp = bpID[ID] + N_BP*strandID[ID] 
     
        #look back until a compatible comparison seq is found
        #i.e. one that does not contain a conflicting segment
        b = i-1
        match_rev = -1
        misprime_rev = -1
        while (b>=0):
            #check if current comp sequence contains incompatable segments   
            #incompatible segments have IDs of same magnitude but oppoisite sign
            comp_seq = segmentID[sort_idx[b]]        
            comp_seq = [x * -1 for x in comp_seq]    
            #if no conflicting segments start comparison. Break cycle once comparison
            #is complete
            if set(segments).isdisjoint(set(comp_seq)):
                count = 0
                misprime_score = 0
                str_1 = seq
                str_2 = subset_str[b]
                #run through strings to see how many consecutive bp's match
                while (count < len(str_1) - 1 and count < len(str_2) - 1): 
                    if str_1[count + 1] != str_2[count + 1]:               
                        break
                    count += 1
        
                    if str_1[count] == 'G' or str_1[count] == 'C':
                        misprime_score += 2.0 #NL: Studies indicate that energy associated with GC bonds is about twice that of AT  
                    else:
                        misprime_score += 1.0
        
                match_rev = count
                misprime_rev = misprime_score
                
                break
            else:
                #if no match, look back to the next one
                b += -1
    
        #look forward
        f = i+1
        match_fwd = -1
        misprime_fwd = -1
        while (f<len(subset_str)):
            #check if current comp sequence contains incompatable segments        
            comp_seq = segmentID[sort_idx[f]]        
            comp_seq = [x * -1 for x in comp_seq]        
            if set(segments).isdisjoint(set(comp_seq)):
                count = 0
                misprime_score = 0
                str_1 = seq
                str_2 = subset_str[f]
        
                while (count < len(str_1) - 1 and count < len(str_2) - 1): 
                    if str_1[count + 1] != str_2[count + 1]:               
                        break
                    count += 1
        
                    if str_1[count] == 'G' or str_1[count] == 'C':
                        misprime_score += 2.0 #NL: Studies indicate that energy associated with GC bonds is about twice that of AT  
                    else:
                        misprime_score += 1.0
        
                match_fwd = count
                misprime_fwd = misprime_score
                break
            else:
                f += 1
                
        if (match_rev > match_fwd):
            match = bpID[sort_idx[b]]
            match_max = match_rev
            misprime_score_max = misprime_rev
        else:
            match = bpID[sort_idx[f]]
            match_max = match_fwd
            misprime_score_max = misprime_fwd
        
        #compare count for current seq to count for position.
        #if higher, take current score
        if (match_max > num_match_mat[0,bp]):
            num_match_mat[0,bp] = match_max
            misprime_score_mat[0,bp] = misprime_score_max
            best_match_mat[0,bp] = match
    
    #split arrays into forward and reverse to match original version output    
    num_match_forward = num_match_mat[:,0:N_BP]
    misprime_score_forward = misprime_score_mat[:,0:N_BP]
    best_match_forward = best_match_mat[:,0:N_BP]
    
    num_match_reverse = num_match_mat[:,N_BP:]
    misprime_score_reverse = misprime_score_mat[:,N_BP:]
    best_match_reverse = best_match_mat[:,N_BP:]
    
    return (num_match_forward, num_match_reverse, best_match_forward, best_match_reverse, misprime_score_forward, misprime_score_reverse, subset_str, sort_idx, bpID, segmentID)


"""
"NL Addition: Adaptation of original misprime function to assess misprime "
             "for chosen set of primers"
"""

def _check_misprime_post(sequence, primer_set, misprime_warn, primer_key):

    N = len(sequence)
    # length of string subsets
    m = 20
    subset_str = []

    # match to sequence
    for i in xrange(N):
        start_pos = max(i - m, 0) - 1
        if start_pos == -1:
            subset_str.append(sequence[i::-1])
        else:
            subset_str.append(sequence[i:start_pos:-1])

    # match to reverse complement of sequence
    sequence_rc = reverse_complement(sequence)
    for i in xrange(N):
        end_pos = N - i - 1
        start_pos = max(end_pos - m - 1, 0) - 1
        if start_pos == -1:
            subset_str.append(sequence_rc[end_pos::-1])
        else:
            subset_str.append(sequence_rc[end_pos:start_pos:-1])

    
    match_list = []
    misprime_score_list = []
    match_id_list = []
    strand = -1
    for i in xrange(len(primer_set)):
        
        str_1 = primer_set[i]
        str_1 = str_1[::-1]
        strand = -1*strand
        iterlist = range(2*N)
        if (strand==1):
            iterlist.remove(primer_key[1,i])
        else:
            iterlist.remove(primer_key[0,i]+N)
        match = []
        misprime = []
        match_id = []
        for j in iterlist :
            ind = j % N
            
            
            count = -1
            misprime_score = 0
            
            str_2 = subset_str[j]
    
            while (count < len(str_1) - 1 and count < len(str_2) - 1): 
                if str_1[count + 1] != str_2[count + 1]:               
                    break
                count += 1
    
                if str_1[count] == 'G' or str_1[count] == 'C':
                    misprime_score += 2.0 #NL: Studies indicate that energy associated with GC bonds is about twice that of AT  
                else:
                    misprime_score += 1.0
    
            if (count > misprime_warn):
                sub_match = []
                for k in xrange(len(primer_key[1,])):
                    if (ind >= primer_key[0,k] and ind <= primer_key[1,k]):
                        sub_match.append(k)   
                        
                match_id.append(sub_match)
                match.append(count)
                misprime.append(misprime_score)                               
                

        match_list.append(match)
        misprime_score_list.append(misprime)
        match_id_list.append(match_id)
                
    return (match_list, misprime_score_list, match_id_list)





