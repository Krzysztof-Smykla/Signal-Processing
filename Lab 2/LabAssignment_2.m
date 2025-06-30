% Assignment #2

% 1. Convolution

%signals x[k] and the impulse response h[k] below use Matlab to obtain the system response y[k]:

% Defining a unit step function u 
% Start with all zeros: 
k = -10:10;
u = @(k) double(k >= 0);

% Point a)
x = u(k) - u(k - 8);
h = sin(2*pi*k/8) .* (u(k) - u(k - 8));

figure(1);
% Input Signal
subplot(2,1,1);
stem(k, x, 'filled', 'LineWidth', 2);
title('x[k]'); xlabel('k'); ylabel('x[k]');
% Impulse Response
subplot(2,1,2);
stem(k, h, 'filled', 'LineWidth', 2);
title('h[k]'); xlabel('k'); ylabel('h[k]');

% Convolution of x(k) * h(k) using modified conv_m() function
y = conv_m(x,k,h,k)

figure(2);
plot(y,"-b","LineWidth",2)
title('y = x[k] * h[k]'); xlabel('k'); ylabel('y[k]')

% Point b)
x2 = sin(2*pi*k/8) .* (u(k) - u(k - 8));
h2 = -sin(2*pi*k/8) .* (u(k) - u(k - 8));

figure(3);
% Input Signal
subplot(2,1,1);
stem(k, x2, 'filled', 'LineWidth', 2);
title('x2[k]'); xlabel('k'); ylabel('x2[k]');
% Impulse Response
subplot(2,1,2);
stem(k, h2, 'filled', 'LineWidth', 2);
title('h2[k]'); xlabel('k'); ylabel('h2[k]');
% Convolution of x2(k) * h2(k)using modified conv_m() function
y2 = conv_m(x2,k,h2,k)

figure(4);
plot(y2,"-b","LineWidth",2)
title('y2 = x2[k] * h2[k]'); xlabel('k'); ylabel('y2[k]')

%%
% 2. Sampling

clear
clearvars

% Signal a) 
Ns = 5000;
%k = 1:Ns;
k = 0:Ns-1;
% Signal x
x = 4*cos(2*pi*0.1*k) +2*cos(2*pi*0.35*k);

figure(1)
% The fourier transform of x
h = fft(x);
f= linspace(0,1000,Ns);
subplot(2,1,1)
% Plotting the signel between 0 and the Nyquist frequency Ns
stem(f,abs(h(1:Ns))/Ns), grid
title('fs = 1000Hz')
xlabel('Frequency in Hz')
ylabel('DFT Amplitude')

subplot(2,1,2)
stem(f,real(h)/Ns)
hold on
stem(f,imag(h)/Ns)
hold off
title('fs = 1000Hz')
xlabel('Frequency in Hz')
ylabel('Amplitude')
legend('Real H','Imag H')
grid
%%
% Decimation of signal x
xd = x(1:2:end); % Decimate signal
Nsd = length(xd);
hd = fft(xd);
fd = linspace(0,500,Nsd);
% Plotting decimated signal
figure(2)
subplot(211)
stem(fd,abs(hd(1:Nsd))/Nsd), grid
title('fs = 500Hz')
xlabel('Frequency in Hz')
ylabel('DFT Amplitude')
subplot(212)
stem(fd,real(hd)/Nsd)
hold on
stem(fd,imag(hd)/Nsd)
hold off
title('fs = 500Hz')
xlabel('Frequency in Hz')
ylabel('Aamplitude')
legend('Real H','Imag H')
grid
%%
% Signal b)
Ns = 5000;
k = 1:Ns;
%k = 0:Ns-1;
x = 4*cos(2*pi*0.1*k) +2*cos(2*pi*0.4*k);

figure(1)
h = fft(x);
f= linspace(0,1000,Ns);
subplot(211)
stem(f,abs(h(1:Ns))/Ns), grid
title('fs = 1000Hz')
xlabel('Frequency in Hz')
ylabel('DFT Amplitude')
subplot(212)
stem(f,real(h)/Ns)
hold on
stem(f,imag(h)/Ns)
hold off
title('fs = 1000Hz')
xlabel('Frequency in Hz')
ylabel('Amplitude')
legend('Real H','Imag H')
grid

% Decimation
xd = x(1:2:end); % Decimate signal
Nsd = length(xd);
hd = fft(xd);
fd = linspace(0,500,Nsd);

figure(2)
subplot(211)
stem(fd,abs(hd(1:Nsd))/Nsd), grid
title('fs = 500Hz')
xlabel('Frequency in Hz')
ylabel('DFT Amplitude')
subplot(212)
stem(fd,real(hd)/Nsd)
hold on
stem(fd,imag(hd)/Nsd)
hold off
title('fs = 500Hz')
xlabel('Frequency in Hz')
ylabel('Aamplitude')
legend('Real H','Imag H')
grid