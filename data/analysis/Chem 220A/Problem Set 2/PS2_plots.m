%plot for 5b
%define alpha array
alpha = linspace(4,20);

energy = arrayfun(@(x) x - coth(1/x), alpha);

plot(alpha, energy)

title('<u>/\epsilon As a Function of kT/\epsilon');

xlabel('kT/\epsilon') % x-axis label
ylabel('<u>/\epsilon') % y-axis label
%%
%plot for e
%define alpha array
alpha1 = linspace(0,100,1000);

energy1 = arrayfun(@(x) x - coth(1/x), alpha1);

plot(alpha1, energy1)

title('<u>/\epsilon As a Function of kT/\epsilon (Limiting Behaviors)');

xlabel('kT/\epsilon') % x-axis label
ylabel('<u>/\epsilon') % y-axis label