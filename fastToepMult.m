function [yf, t] = fastToepMult(a, b, w)
%FASTTOEPMULT O(nlog(n))y = T*w multiplication
%   Arguments:
%   a:  the first vector of the toeplitz T matrix
%   b:  the second vector of the toeplitz T matrix
%   w:  a vector
%   Output:
%   yf: result
%   t:  calculation time (sec)
%%

% Start timer
tic

a = a(:);
b = b(:)';
b = b(2:end);

len_a = length(a);
c = [a; 0; fliplr(b)'];

p = ifft(fft(c).*fft([w; zeros(len_a, 1)]));
yf = p(1:len_a);

% Return time
t = toc;

end

