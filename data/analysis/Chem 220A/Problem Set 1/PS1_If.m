%%%plot large deviation form of state probability (from PS1 #1.d)
function [prob] = PS1_If(f)
    prob = f*log(f) + (1-f)*log(1-f) + log(2);
end



