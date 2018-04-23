function [t, tf, er] = part1_comp(u, w)
%PART1_COMP helper script for the question 1 
%% --- Variable set up ---
n = length(u);
m = length(w);
normal_calculation_limit = 8001;

%% --- Calculate product normally, time the calculation ---
% Else it takes too much memory for the T vector
if n < normal_calculation_limit
    T = toeplitz( u(m:n), u(m:-1:1) );
    tic
    y = T*w;
    t = toc;
end

%% --- Calculate the product with the O(nlogn) method ---
[yf, tf] = fastToepMult( u(m:n), u(m:-1:1), w);

%% --- Error ---
if n < normal_calculation_limit;
    er = norm(yf-y)/norm(y);
end

end
