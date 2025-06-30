
% Plot the attenuation in decibels vs relative frequency f/fo.
% Compare responses of the filters and comment differences.
clear
clc

% Filter parameters
n = 3;                     % Filter order
fo = 1e5;                  % Cutoff frequency in Hz
wc = 2 * pi * fo;          % Angular frequency in rad/s
w = logspace(3, 6, 4096);  % Frequency vector: 1 kHz to 1 MHz

%% 1. Butterworth Filter
[zb, pb, kb] = butter(n, wc, 's');
[bb, ab] = zp2tf(zb, pb, kb);
[hb, ~] = freqs(bb, ab, w);

%% 2. Chebyshev Type I Filter (3 dB ripple)
[z1, p1, k1] = cheby1(n, 3, wc, 's');
[b1, a1] = zp2tf(z1, p1, k1);
[h1, ~] = freqs(b1, a1, w);

%% 3. Chebyshev Type II Filter (30 dB stopband attenuation)
[z2, p2, k2] = cheby2(n, 30, wc, 's');
[b2, a2] = zp2tf(z2, p2, k2);
[h2, ~] = freqs(b2, a2, w);

%% 4. Elliptic Filter (3 dB ripple, 30 dB stopband attenuation)
[ze, pe, ke] = ellip(n, 3, 30, wc, 's');
[be, ae] = zp2tf(ze, pe, ke);
[he, ~] = freqs(be, ae, w);

%% Plot Magnitude Responses
% Reference angular frequency
wo = 2 * pi * fo;

% Plot 1: Attenuation (dB) vs Relative Frequency (f/fo)
figure(1)
plot(w / wo, mag2db(abs(hb)), 'b', 'LineWidth', 1.5); hold on
plot(w / wo, mag2db(abs(h1)), 'r--', 'LineWidth', 1.5);
plot(w / wo, mag2db(abs(h2)), '-', 'Color', [0.5 0 0.5], 'LineWidth', 1.5);
plot(w / wo, mag2db(abs(he)), 'g-.', 'LineWidth', 1.5);
axis([0 2 -80 5])
grid on
title('Filter Gain', 'FontWeight', 'bold')
xlabel('Relative Frequency f/f_o', 'FontWeight', 'bold')
ylabel('Attenuation (dB)', 'FontWeight', 'bold')
legend('Butterworth','Chebyshev I','Chebyshev II','Elliptic', 'Location', 'SouthWest')
set(gca, 'FontSize', 12)

% Plot 2: Phase vs Relative Frequency (normalized by π)
figure(2)
plot(w / wo, unwrap(angle(hb)) / pi, 'b', 'LineWidth', 1.5); hold on
plot(w / wo, unwrap(angle(h1)) / pi, 'r--', 'LineWidth', 1.5);
plot(w / wo, unwrap(angle(h2)) / pi, '-', 'Color', [0.5 0 0.5], 'LineWidth', 1.5);
plot(w / wo, unwrap(angle(he)) / pi, 'g-.', 'LineWidth', 1.5);
axis([0 2 -1.1 1.1])
grid on
title('Filter Phase', 'FontWeight', 'bold')
xlabel('Relative Frequency f/f_o', 'FontWeight', 'bold')
ylabel('Phase (×π radians)', 'FontWeight', 'bold')
legend('Butterworth','Chebyshev I','Chebyshev II','Elliptic', 'Location', 'SouthWest')
set(gca, 'FontSize', 12)