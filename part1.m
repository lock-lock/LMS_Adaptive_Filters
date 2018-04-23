%% --- What is this script ---
% This script compares the multiplication time of: y = T*w, where
%   T is a toeplitz matrix
%   w is a vector
% The multiplication is calculated using the usual method and the fast
% toeplitz multiplication method (all this happen inside part_1_comp
% function.

%% --- Calculate times and error ---
normal_t = zeros(4000, 1);
fast_t = zeros(4000, 1);
error = zeros(4000, 1);

parfor i=1:4000
    n = i*2;
    m = i;
    u = randn(n, 1);
    w = randn(m, 1);
    
    [t, tf, er] = part1_comp(u, w);
    normal_t(i) = t;
    fast_t(i) = tf;
    error(i) = er;
end

%% --- Plot time results ---
n = 2:2:8000;
figure('name', 'Time comparison')
hold on
plot(n, fast_t);
plot(n, normal_t);
title({'y = T*w: Time comparison of the normal and the O(nlogn) algorithms',...
    '{T comes from an n-vector and w is a n/2-vector}'})
ylabel('Time (s)')
xlabel('n')
legend('O(nlogn) algorithm', 'common algorithm')

%% --- Plot error ---
figure('name', 'Error')
plot(n, error);
title({'y = T*w: Error of the O(nlogn) algorithm',...
    '{T comes from an n-vector and w is a n/2-vector}'})
ylabel('Error')
xlabel('n')


%% 
%   Author: Nikolaos Katomeris, 8551, ngkatomer@auth.gr
%   Last change at: April 21st, 2018