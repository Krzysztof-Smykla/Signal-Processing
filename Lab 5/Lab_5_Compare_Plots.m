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
figure
semilogx(w/(2*pi), 20*log10(abs(hb)), 'b-', 'LineWidth', 2); hold on
semilogx(w/(2*pi), 20*log10(abs(h1)), 'r--', 'LineWidth', 2);
semilogx(w/(2*pi), 20*log10(abs(h2)), '-', 'Color', [0.5 0 0.5], 'LineWidth', 2);
semilogx(w/(2*pi), 20*log10(abs(he)), 'g-.', 'LineWidth', 2);
grid on
xlabel('Frequency (Hz)', 'FontWeight', 'bold')
ylabel('Magnitude (dB)', 'FontWeight', 'bold')
title('Magnitude Response of Analog Lowpass Filters', 'FontWeight', 'bold')
legend('Butterworth', 'Chebyshev Type I', 'Chebyshev Type II', 'Elliptic', 'Location', 'SouthWest')
set(gca, 'FontSize', 12)

%% Plot Phase Responses
figure
semilogx(w/(2*pi), unwrap(angle(hb)) * 180/pi, 'b-', 'LineWidth', 2); hold on
semilogx(w/(2*pi), unwrap(angle(h1)) * 180/pi, 'r--', 'LineWidth', 2);
semilogx(w/(2*pi), unwrap(angle(h2)) * 180/pi, '-', 'Color', [0.5 0 0.5], 'LineWidth', 2);
semilogx(w/(2*pi), unwrap(angle(he)) * 180/pi, 'g-.', 'LineWidth', 2);
grid on
xlabel('Frequency (Hz)', 'FontWeight', 'bold')
ylabel('Phase (degrees)', 'FontWeight', 'bold')
title('Phase Response of Analog Lowpass Filters', 'FontWeight', 'bold')
legend('Butterworth', 'Chebyshev Type I', 'Chebyshev Type II', 'Elliptic', 'Location', 'SouthWest')
set(gca, 'FontSize', 12)