%% --- What is this script ---
% This script proves the convolution theorem
% It calculates the convolution of 2 vectors in different ways
% And shows that they are all equivalent as convolution theorem states.

%% --- Prepare variables ---
N = 500;
x = rand(N, 1);
y = rand(N, 1);

%% --- Calculate convolution with conv() ---
r1 = conv(x, y);

%% --- Calculate convolution using a Toeplitz matrix ---
% r2 = Y*x, where Y is toeplitz of y
zy = padarray(y, N-1, 0, 'post');
cv = [zy(1) zeros(1, N-1)];
Y = toeplitz(zy, cv);
r2 = Y*x;

%% --- Calculate convolution using a Circulant matrix ---
% r3 = C*x, where C is a circulant extension of the Y matrix
zx = padarray(x, N-1, 0, 'post');
C = [Y, toeplitz(circshift(zy, -N + 1), zy(N:-1:2))];
r3 = C*zx;

%% --- Calculate convolution in the frequency domain ---
r4 = ifft(fft(zy).*fft(zx));

%% --- Plot calculation differences ---
figure('name', 'convolutions')
hold on
plot(r4-r1)
plot(r3-r1)
plot(r2-r1)
legend( 'y_4 = ifft(fft(y)*fft(x))', 'y_3 = Cx', 'y_2 = Yx', 'location', 'best')
title({'Difference of the convolution result', 'y_1 = conv(x, y)'})
ylabel('y_1 - y_i')

%% --- File info ---
%   Author: Nikolaos Katomeris, 8551, ngkatomer@auth.gr
%   Last change at: April 23th, 2018