%
% Design a nth-order analog Butterworth lowpass filter
% with a cutoff frequency of fo Hz.
% Multiply by 2pi to convert the frequency to radians per second.
% Compute the frequency response of the filter at 4096 points.
%
%%
clear
clc

% Set filter order and cutoff frequency
n = 3;                 % Filter order
fo = 1e5;              % Cutoff frequency in Hz (100,000 Hz)
wc = 2 * pi * fo;      % Angular frequency in rad/s

% Design an analog Butterworth low-pass filter
[zb, pb, kb] = butter(n, wc, 's');  % 's' specifies analog filter
[bb, ab] = zp2tf(zb, pb, kb);       % Convert from ZPK to transfer function

% Compute frequency response
w = logspace(3, 6, 4096);           % Frequency range from 1 kHz to 1 MHz
[hb, wb] = freqs(bb, ab, w);

figure(1)
semilogx(wb/(2*pi), 20*log10(abs(hb)));
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Butterworth Lowpass Filter (n = 3)')

figure(2)
semilogx(wb/(2*pi), unwrap(angle(hb)) * 180/pi, Color='red');
grid on
xlabel('Frequency (Hz)')
ylabel('Phase (degrees)')
title('Butterworth Lowpass Filter phase')

%%
% Set filter order and cutoff frequency
m = 7;                 % Filter order
fo = 1e5;              % Cutoff frequency in Hz (100,000 Hz)
wc = 2 * pi * fo;      % Angular frequency in rad/s

% Design an analog Butterworth low-pass filter
[zb, pb, kb] = butter(m, wc, 's');  % 's' specifies analog filter
[bb, ab] = zp2tf(zb, pb, kb);       % Convert from ZPK to transfer function

% Compute frequency response
w = logspace(3, 6, 4096);           % Frequency range from 1 kHz to 1 MHz
[hb, wb] = freqs(bb, ab, w);

figure(3)
semilogx(wb/(2*pi), 20*log10(abs(hb)));
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Butterworth Lowpass Filter (n = 7)')

figure(4)
semilogx(wb/(2*pi), unwrap(angle(hb)) * 180/pi, Color='red');
grid on
xlabel('Frequency (Hz)')
ylabel('Phase (degrees)')
title('Butterworth Lowpass Filter phase')


%%
% Design a nth-order Chebyshev Type I filter with the same edge frequency
% and 3 dB of passband ripple.
% n = 3

[z1,p1,k1] = cheby1(n,3,2*pi*fo,"s");
[bb, ab] = zp2tf(z1, p1, k1);  % <-- WRONG

% Compute frequency response
w = logspace(3, 6, 4096);           % Frequency range from 1 kHz to 1 MHz
[hb, wb] = freqs(bb, ab, w);

figure(5)
semilogx(wb/(2*pi), 20*log10(abs(hb)))
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Cheb1 Lowpass Filter (n = 3)')

figure(6)
semilogx(wb/(2*pi), unwrap(angle(hb)) * 180/pi, Color='red');
grid on
xlabel('Frequency (Hz)')
ylabel('Phase (degrees)')
title('Cheb1 Lowpass Filter phase')
%%
% Design a nth-order Chebyshev Type I filter with the same edge frequency
% and 3 dB of passband ripple.
% Filter order  n = 7
m = 7;  
[z1,p1,k1] = cheby1(m,3,2*pi*fo,"s");
[bb, ab] = zp2tf(z1, p1, k1);       % Convert from ZPK to transfer function

% Compute frequency response
w = logspace(3, 6, 4096);           % Frequency range from 1 kHz to 1 MHz
[hb, wb] = freqs(bb, ab, w);

figure(7)
semilogx(wb/(2*pi), 20*log10(abs(hb)))
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Cheb1 Lowpass Filter (n = 7)')

figure(8)
semilogx(wb/(2*pi), unwrap(angle(hb)) * 180/pi, Color='red');
grid on
xlabel('Frequency (Hz)')
ylabel('Phase (degrees)')
title('Cheb1 Lowpass Filter phase')

%%
% Design a nth-order Chebyshev Type II filter with the same edge frequency
% and 3 dB of passband ripple.
% n = 3
[z2,p2,k2] = cheby2(n,30,2*pi*fo,"s");
[bb, ab] = zp2tf(z2, p2, k2);       % Convert from ZPK to transfer function

% Compute frequency response
w = logspace(3, 6, 4096);           % Frequency range from 1 kHz to 1 MHz
[hb, wb] = freqs(bb, ab, w);

figure(9)
semilogx(wb/(2*pi), 20*log10(abs(hb)))
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Cheb2 Lowpass Filter (n = 3)')

figure(10)
semilogx(wb/(2*pi), unwrap(angle(hb)) * 180/pi, Color='red');
grid on
xlabel('Frequency (Hz)')
ylabel('Phase (degrees)')
title('Cheb2 Lowpass Filter phase')
%%
% Design a nth-order Chebyshev Type II filter with the same edge frequency
% and 3 dB of passband ripple.
% Filter order  n = 7
m = 7;  
[z2,p2,k2] = cheby2(m,30,2*pi*fo,"s");
[bb, ab] = zp2tf(z2, p2, k2);       % Convert from ZPK to transfer function

% Compute frequency response
w = logspace(3, 6, 4096);           % Frequency range from 1 kHz to 1 MHz
[hb, wb] = freqs(bb, ab, w);

figure(11)
semilogx(wb/(2*pi), 20*log10(abs(hb)))
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Cheb2 Lowpass Filter (n = 7)')

figure(12)
semilogx(wb/(2*pi), unwrap(angle(hb)) * 180/pi, Color='red');
grid on
xlabel('Frequency (Hz)')
ylabel('Phase (degrees)')
title('Cheb2 Lowpass Filter phase')
%%
% Design a nth-order elliptic filter with the same edge frequency, 3 dB of passband ripple,
% and 30 dB of stopband attenuation. Compute its frequency response.
% n = 3
[ze,pe,ke] = ellip(n,3,30,2*pi*fo,"s");
% Compute its frequency response.
[bb, ab] = zp2tf(ze, pe, ke);       % Convert from ZPK to transfer function

w = logspace(3, 6, 4096);           % Frequency range from 1 kHz to 1 MHz
[hb, wb] = freqs(bb, ab, w);

figure(13)
semilogx(wb/(2*pi), 20*log10(abs(hb)))
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Eliptic Lowpass Filter (n = 3)')

figure(14)
semilogx(wb/(2*pi), unwrap(angle(hb)) * 180/pi, Color='red');
grid on
xlabel('Frequency (Hz)')
ylabel('Phase (degrees)')
title('Eliptic Lowpass Filter phase')

%%
% Design a nth-order elliptic filter with the same edge frequency, 3 dB of passband ripple,
% and 30 dB of stopband attenuation. Compute its frequency response.
% n = 7
[ze,pe,ke] = ellip(m,3,30,2*pi*fo,"s");
% Compute its frequency response.
[bb, ab] = zp2tf(ze, pe, ke);       % Convert from ZPK to transfer function

w = logspace(3, 6, 4096);           % Frequency range from 1 kHz to 1 MHz
[hb, wb] = freqs(bb, ab, w);

figure(15)
semilogx(wb/(2*pi), 20*log10(abs(hb)))
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Eliptic Lowpass Filter (n = 7)')

figure(16)
semilogx(wb/(2*pi), unwrap(angle(hb)) * 180/pi, Color='red');
grid on
xlabel('Frequency (Hz)')
ylabel('Phase (degrees)')
title('Eliptic Lowpass Filter phase')
%%
% Design a nth-order Bessel filter with the same edge frequency, 3 dB of passband ripple,
% and 30 dB of stopband attenuation. Compute its frequency response.
n = 3 
fo = 1e5;              % Cutoff frequency in Hz (100,000 Hz)
wc = 2 * pi * fo;      % Angular frequency in rad/s

[bb, ab] = besself(n,fo);
% Compute its frequency response.

w = logspace(3, 6, 4096);           % Frequency range from 1 kHz to 1 MHz
[hb, wb] = freqs(bb, ab, w);

figure(17)
semilogx(wb/(2*pi), 20*log10(abs(hb)))
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Bessel Lowpass Filter (n = 3)')

figure(18)
semilogx(wb/(2*pi), unwrap(angle(hb)) * 180/pi, Color='red');
grid on
xlabel('Frequency (Hz)')
ylabel('Phase (degrees)')
title('Bessel Lowpass Filter phase')

%%
% Design a nth-order Bessel filter with the same edge frequency, 3 dB of passband ripple,
% and 30 dB of stopband attenuation. Compute its frequency response.
m = 7
[bb, ab] = besself(m,fo);

w = logspace(3, 6, 4096);           % Frequency range from 1 kHz to 1 MHz
[hb, wb] = freqs(bb, ab, w);

figure(19)
semilogx(wb/(2*pi), 20*log10(abs(hb)))
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Bessel Lowpass Filter (n = 7)')

figure(20)
semilogx(wb/(2*pi), unwrap(angle(hb)) * 180/pi, Color='red');
grid on
xlabel('Frequency (Hz)')
ylabel('Phase (degrees)')
title('Bessel Lowpass Filter phase')