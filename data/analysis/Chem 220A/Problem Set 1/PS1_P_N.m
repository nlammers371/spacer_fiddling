%plot P(N) for varying valued of M
function prob = PS1_P_N(f,M)
    prob = exp(-M*PS1_I_f(f));
end