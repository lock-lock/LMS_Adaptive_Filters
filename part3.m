%% --- What is this script ---
% This script uses an adaptive filter in order to model an unknown linear
% system. It does that in 4 differect ways

%% --- Variable preparation ---
mu_default = 0.005;
M = 2^10;
N = 200*M;
% v: white noise with var = 0.57
% u = v
u = sqrt(0.57)*rand(N+M, 1); u = u - mean(u);

% d is the unknown system's output
d = plant(u')';

%% --- Block LMS algorithm with two loops ---
tic
% Variable preparation
mu = mu_default;    % mu: adjastive parameter
L = M;              % L: block size (same as M, see pg 10 of 6th lecture)
w = zeros(N/L, M);  % w: filter parameters
y = zeros(N, 1);    % y: filter output
e = d(M+1:end);     % e: error

% Algorithm
for k = 1:N/L-1
    % initialize phi
    phi = zeros(M, 1);
    
    for i = 1:L
        y(k*L+i) = w(k, :)*u(M+k*L+i:-1:k*L+i+1); % filter output
        e(k*L+i) = d(M+k*L+i) - y(k*L+i);
        phi = phi+mu*e(k*L+i)*u(M+k*L+i:-1:k*L+i+1);
    end
    w(k+1, :) = w(k, :) + phi';
end
e1 = e;
w1 = w(end, :);
t1 = toc;
%% --- Block LMS algorithm with one loop ---
tic
% Variable preparation
mu = mu_default;    % mu: adjastive parameter
L = M;              % L: block size
w = zeros(N/L, M);  % w: filter parameters
e = d(M+1:end);     % e: error

% Algorithm
for k = 1:N/L-1
    
    % All the inputs needed for this block
    u_b = toeplitz(u(k*L:1:(k+1)*L-1),u(k*L:-1:(k-1)*L+1));
    % Desired output
    d_b = d(k*L:(k+1)*L-1);
    % Filter output
    y = u_b*w(k, :)';
    % Error
    e(k*L+1:(k+1)*L) = d_b - y;
    phi = (u_b.'*e(k*L+1:(k+1)*L))';
    % New filter parameters
    w(k+1, :) = w(k, :) + mu*phi;
end
e2 = e;
w2 = w(end, :);
t2 = toc;
%% --- Block LMS algorithm using (frequency domain - fast LMS) ---
tic
% Variable preparation
mu = mu_default;    % mu: adjastive parameter
L = M;              % L: block size
w = zeros(M, 1);    % w: filter parameters
e = d(M+1:end);     % e: error

for k = 1:N/L-1
    % Fourier of the extended filter parameters
    W = fft([w; zeros(L,1)]);
    % Fourier of the input
    U = fft(u((k-1)*L+1:(k+1)*M), 2*M);
    % Filter output
    y_b = ifft(U.*W);
    y_b = y_b(L+1:2*L);
    % Desired output and error
    d_b = d(k*L+1:(k+1)*L);
    e (k*L+1:(k+1)*L) = d_b - y_b; 
    
    % e is padded at the beginning with L zeros and then a 2*L-length FFT
    % is computed
    E = fft([zeros(L, 1); e(k*L+1:(k+1)*L)], 2*L);
    
    % helper = [phi, D];
    helper = ifft(E.*conj(U));
    % helper2 = fft of helper withoud D.
    helper2 = fft([helper(1:M); zeros(M,1)]);
    
    % New filter parameters
    W = W + mu*helper2;
    w = ifft(W);
    w = w(1:M);
end
e3 = e;
w3 = ifft(W);
w3 = real(w3(1:M))';
t3 = toc;

%% --- Block LMS algorithm using (frequency domain - unconstrained) ---
tic
% Variable preparation
mu = mu_default;    % mu: adjastive parameter
L = M;              % L: block size
e = d(M+1:end);     % e: error

% Fourier of the extended filter parameters: initialization
W = zeros(2*L, 1);
for k = 1:N/L-1
    % Fourier of the input
    U = fft(u((k-1)*L+1:(k+1)*M), 2*M);
    % Filter output
    y_b = ifft(U.*W);
    y_b = y_b(L+1:2*L);
    % Desired output and error
    d_b = d(k*L+1:(k+1)*L);
    e (k*L+1:(k+1)*L) = d_b - y_b; 
    
    % e is padded at the beginning with L zeros and then a 2*L-length FFT
    % is computed
    E = fft([zeros(L, 1); e(k*L+1:(k+1)*L)], 2*L);
    
    % helper = [phi, D];
    helper = conj(U)'.*E';
    
    % New filter parameters
    W = W + mu*helper';
end
e4 = e;
w4 = ifft(W);
w4 = real(w4(1:M))';
t4 = toc;

%% --- Clean up ---
clear helper helper2 i k L M phi w W e R y_b d_b y U u_b E

%% --- Plot comparison plots ---
% Times
times = [t1 t2 t3 t4];
labels = categorical({'1st', '2nd', '3rd', '4th'});
figure('name', 'time comparison')
bar(labels, times);
xlabel('algorithm')
ylabel('time (s)')
title('Time comparison of the 4 algorithms')

% Adaptation
n = 1:N;
figure('name', 'adaptation plots')
subplot(2,2,1);
plot(n, e1.^1)
xlabel('step')
ylabel('Output Error')
title('1st algorithm')
xlim([0 N])

subplot(2,2,2);
plot(n, e3.^1)
xlabel('step')
ylabel('Output Error')
title('2nd algorithm')
xlim([0 N])

subplot(2,2,3);
plot(n, e3.^1)
xlabel('step')
ylabel('Output Error')
title('3rd algorithm')
xlim([0 N])

subplot(2,2,4);
plot(n, e4.^1)
xlabel('step')
ylabel('Output Error')
title('4th algorithm')
xlim([0 N])

% Learning Curves
figure('name', 'learning curves')
semilogy(n, e1.^2)
hold on
semilogy(n, e2.^2)
semilogy(n, e3.^2)
semilogy(n, e4.^2)
legend('1st algorithm', '2nd algorithm', '3rd algorithm', '4th algorithm', 'location', 'best')
title(['Learning curves for ì = ' num2str(mu)])
ylabel('E(e^2)')
xlabel('step')
xlim([0 N])
%% --- File info ---
%   Author: Nikolaos Katomeris
%   Last change at: April 23th, 2018