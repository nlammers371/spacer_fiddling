%Plot P(N) for various values of M

N = linspace(0,1);

for i = 1:100
    S10(i) = PS1_P_N(N(i), 10);
end
for i = 1:100
    S100(i) = PS1_P_N(N(i), 100);
end
for i = 1:100
    S1000(i) = PS1_P_N(N(i), 1000);
end
for i = 1:100
    S10000(i) = PS1_P_N(N(i), 10000);
end
plot(N,S10,N,S100,N,S1000,N,S10000)

title('PS1 3.g: Plot of P(phi) vs. phi for Selected M Values')
xlabel('phi (N/M)') % x-axis label
ylabel('P(phi)') % y-axis label

legend('M = 10','M = 100','M = 1000','M = 10000','Location','northeast')
