f = linspace(0,1);
for i = 1:100
    F(i) = PS1_I_f(i/100);
end

plot(f,F)
title('PS1 1.d: Plot of I(f) vs. f')
xlabel('f (N/M)') % x-axis label
ylabel('I(f)') % y-axis label