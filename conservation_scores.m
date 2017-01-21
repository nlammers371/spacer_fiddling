%find swap segments that are most highly conserved for initial testing

%parameters
%offset btw full es2 and minimal seq
min_offset = 50;
%length of minimal es2
len_min = 484;
%import conservation score file
cons_file = 'C:\Users\Nicholas\Documents\GitHub\synth_es2\data\analysis\ToFIMO\Results_12.01.16\SegmentedES2_CONSscores.csv';
cons_scores = importdata(cons_file);
cons_scores = cons_scores.data;
min_scores = cons_scores(min_offset:(len_min+min_offset-1),2);
%Import segment indices 
segment_file = 'C:\Users\Nicholas\Documents\GitHub\spacer_fiddling\Primerize-master\results\2017-01-16_11-13-23\segment_index.txt';
segments = importdata(segment_file);
segments = str2num(segments{1});
segment_array = zeros(3,length(segments)/3);
for i = 1:length(segments)
    row = ceil((3*i)/length(segments));
    col = mod(i-1,length(segments)/3) + 1;
    segment_array(row,col) = segments(i);
end
seg_id = zeros(1,len_min);
ct = 1;
for i = 1:size(segment_array,2)
    seg = segment_array(3,i);
    for j = (segment_array(1,i)+1):segment_array(2,i)
        if seg < 0 
            seg_id(ct) = .05;
        else
            seg_id(ct) = .95;
        end
        ct = ct + 1;
    end
end
hold on
plot(1:484, min_scores);
stairs(1:484, seg_id)

axis([0 484 0 1]);
            
            
    