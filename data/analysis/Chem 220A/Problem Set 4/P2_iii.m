f = linspace(.0001,.9999);
N = 100000;

S = (-log(1-f) + f.*log(1./f-1));

plot(f,S)
title('P2 part iii: Normalized Entropy vs. f');
xlabel('f (n/N)');
ylabel('Normalized Entropy (S/(k*N))');