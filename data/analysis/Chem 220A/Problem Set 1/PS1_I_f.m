%%%plot large deviation form of state probability (from PS1 #1.d)
function prob = PS1_I_f(f)
    prob = f*log(f) + (1-f)*log(1-f) + log(.8);
end



