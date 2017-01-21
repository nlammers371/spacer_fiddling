max_L = 60*2^5;
log_max_L = log(max_L);

misprime_warn = 5;
misprime_base = log_max_L/(2*(1.25*(misprime_warn+1)));

l_scores1 = zeros(1,(2*(misprime_warn+1))+1);
l_scores2 = zeros(1,(2*(misprime_warn+1))+1);
l_scores3 = zeros(1,(2*(misprime_warn+1))+1);
prime_scores = zeros(1,(2*(misprime_warn+1))+1);

for i = 0:(2*(misprime_warn+1))
%     L = i*floor(max_L/(2*(misprime_warn+1)));
    
    l_scores1(i+1) = log(15)+log(2^(i/2));
    l_scores2(i+1) = log(45)+log(2^(i/2));
    l_scores3(i+1) = log(60)+log(2^(i/2));
    
    prime_scores(i+1) = 1.25*i*misprime_base;
end

