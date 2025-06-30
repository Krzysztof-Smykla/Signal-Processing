% Lab Assignment 3
% 1. Fourier transform

% Task 1.
% Determine the discrete-time Fourier transform of 
% x(n) = (0.5)^n * u(n).

w = (0:1:500)*pi/500; % [0, pi] axis divided into 501 points.
X = exp(1i*w) ./ (exp(1i*w) - 0.5*ones(1,501));
magX = abs(X); angX = angle(X); realX = real(X); imagX = imag(X);
subplot(2,2,1); plot(w/pi,magX); grid
xlabel('frequency in \pi units'); title('Magnitude Part');
ylabel('Magnitude')
subplot(2,2,3); plot(w/pi,angX); grid
xlabel('frequency in \pi units'); title('Angle Part'); ylabel('Radians')
subplot(2,2,2); plot(w/pi,realX); grid
xlabel('frequency in \pi units'); title('Real Part'); ylabel('Real')
subplot(2,2,4); plot(w/pi,imagX); grid
xlabel('frequency in \pi units'); title('Imaginary Part');
ylabel('Imaginary')

%%
% Task 2.
%a)

clear
n = 0:4; x = 1:5; k = 0:500;
w = (pi/500)*k; % Frequency range 0 to pi
X = x * (exp(-1i*pi/500)) .^ (n'*k);
magX = abs(X); angX = angle(X);
realX = real(X); imagX = imag(X);
figure(1)
subplot(2,2,1); plot(k/500,magX);grid
xlabel('frequency in \pi units'); title('Magnitude Part')
subplot(2,2,3); plot(k/500,angX/pi);grid
xlabel('requency in \pi units'); title('Angle Part')
subplot(2,2,2); plot(k/500,realX);grid
xlabel('frequency in \pi units'); title('Real Part')
subplot(2,2,4); plot(k/500,imagX);grid
xlabel('frequency in \pi units'); title('Imaginary Part')

%%
%b
x = 1:5; k = 0:500; w = (pi/500)*k;
X = fft(x,1000); % DFT calculated with zero-padding from O to 2pi
X=X(1:501); % Plot half spectrum 0 to pi
magX = abs(X); angX = angle(X);
realX = real(X); imagX = imag(X);
figure(2)
subplot(2,2,1); plot(k/500,magX);grid
xlabel('frequency in \pi units'); title('Magnitude Part')
subplot(2,2,3); plot(k/500,angX/pi);grid
xlabel('requency in \pi units'); title('Angle Part')
subplot(2,2,2); plot(k/500,realX);grid
xlabel('frequency in \pi units'); title('Real Part')
subplot(2,2,4); plot(k/500,imagX);grid
xlabel('frequency in \pi units'); title('Imaginary Part')

% The plots are in parts a) and b) are identical. 
% Part a calculates the DFT manually while part b uses the built-in fft()
% function in MATLAB.


%%
% Task 3

w = (0:1:500)*pi/500; % [0, pi] axis divided into 501 points.
H = exp(1i*w) ./ (exp(1i*w) - 0.9*ones(1,501));
magH = abs(H); angH = angle(H);
subplot(2,1,1); plot(w/pi,magH); grid;
xlabel('frequency in pi units'); ylabel('|H|');
title('Magnitude Response');
subplot(2,1,2); plot(w/pi,angH/pi); grid
xlabel('frequency in pi units'); ylabel('Phase in pi Radians');
title('Phase Response');
