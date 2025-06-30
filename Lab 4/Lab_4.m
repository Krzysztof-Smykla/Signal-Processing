%% Assignment 4.1
% Task 1: Calculate and plot h[n] using Fourier series method
N = 15;                  % Filter length
n = 0:N-1;               % Sample indices
Wc = 0.3 * pi;           % Cutoff frequency

% Ideal impulse response (centered around (N-1)/2 for linear phase)
n0 = (N - 1) / 2;
h = (Wc / pi) * sinc((Wc / pi) * (n - n0));  % h[n] from Fourier series

% Plot h[n]
figure(1);
stem(n, h, 'filled');
xlabel('n'); ylabel('h[n]');
title('FIR Lowpass Filter Impulse Response h[n]');
grid on;

% Task 2: Calculate H[W] using FFT
N_fft = 512;                
H = fft(h, N_fft);           % Frequency response
H_shifted = fftshift(H);     % Center zero frequency
H_mag = abs(H_shifted);      % Magnitude
H_phase = angle(H_shifted); % Phase

% Frequency axis from -pi to pi
omega = linspace(-pi, pi, N_fft);  

% Task 3: Plot H[W]
figure(2);
plot(omega, H_mag, 'b', 'LineWidth', 0.5);
xlabel('\omega (Rad)'); ylabel('|H[\omega]|');
title('Magnitude Response of FIR Filter');
grid on;

% Phase Response
figure(3);
plot(omega/pi, unwrap(H_phase), 'r');
xlabel('\omega / \pi'); ylabel('Phase of H[\omega]');
title('Phase Response of FIR Filter');
grid on;
%% Assignment 4.2: Window Method FIR Filter Design

N = 15;                      % Total window length
M = (N - 1) / 2;             % Half-length for symmetric window
n = -M:M;                    % Symmetric sample indices

% Hanning window (not Hamming!)
w_hanning = 0.5 - 0.5 * cos(pi * (n + M) / M);

% Blackman window
w_blackman = 0.42 - 0.5 * cos(pi * (n + M) / M) + 0.08 * cos(2 * pi * (n + M) / M);

% Plotting both for comparison with the reference image
figure(1);
stem(n, w_hanning, 'b', 'filled'); grid on;
title('Hamming window w[n]');
xlabel('Index n'); ylabel('Amplitude');
hold on;
figure(2);
stem(n, w_blackman, 'b', 'filled'); grid on;
title('Blackman window w[n]');
xlabel('Index n'); ylabel('Amplitude');

% Print coefficients
disp('Hanning window coefficients w[n]:');
disp(w_hanning);

disp('Blackman window coefficients w[n]:');
disp(w_blackman);
%% Plotting Hanning and Blackman Windows on the Same Plot

N = 15;                         % Window length
n = 0:N-1;                      % Sample indices

% Generate the windows using built-in functions
w_hanning = hann(N);            % Hanning (Hann) window
w_blackman = blackman(N);       % Blackman window

% Plot both windows using line plots
figure;
plot(n, w_hanning, 'b-o', 'LineWidth', 1.5); hold on;
plot(n, w_blackman, 'r-s', 'LineWidth', 1.5);
title('Hanning and Blackman Windows (Length 15)');
xlabel('Index n');
ylabel('Amplitude');
legend('Hanning Window', 'Blackman Window');
grid on;
%% Assignment 4.2 cont

% Calculate and plot h[n] truncted with:
% 1. Rectangular window
% 2. Hanning window
% 3. Blackman window
% Calculate the respective H[W] using fft
% Plot obtained results
% Compare amplitude responses of the designed filers, comment differences

%% FIR Lowpass Filter Design using Window Method

N = 15;                      % Filter length
n = 0:N-1;                   % Sample indices
Wc = 0.3 * pi;               % Cutoff frequency
n0 = (N - 1) / 2;            % Center for symmetry

%% Ideal Impulse Response (Sinc function)
hd = (Wc / pi) * sinc((Wc / pi) * (n - n0));  % Ideal lowpass

%% Window Definitions
w_rect      = ones(1, N);           % Rectangular window
w_hanning   = hann(N)';             % Hanning window
w_blackman  = blackman(N)';         % Blackman window

%% Apply Windows to the Ideal Impulse Response
h_rect     = hd .* w_rect;
h_hanning  = hd .* w_hanning;
h_blackman = hd .* w_blackman;

%% Compute Frequency Response using FFT
N_fft = 512;  % FFT length for better resolution

H_rect     = abs(fft(h_rect, N_fft));
H_hanning  = abs(fft(h_hanning, N_fft));
H_blackman = abs(fft(h_blackman, N_fft));

w = linspace(0, pi, N_fft/2);  % Frequency vector (0 to π)

%% Plot Time-Domain Responses
figure(1);
stem(n, h_rect); title('h[n] with Rectangular Window'); grid on;
figure(2);
stem(n, h_hanning); title('h[n] with Hamming Window'); grid on;
figure(3);
stem(n, h_blackman); title('h[n] with Blackman Window'); grid on;

%% Plot Frequency-Domain (Magnitude) Responses
figure;
plot(w, 20*log10(H_rect(1:N_fft/2)), 'k', 'LineWidth', 1.5); hold on;
plot(w, 20*log10(H_hanning(1:N_fft/2)), 'b', 'LineWidth', 1.5);
plot(w, 20*log10(H_blackman(1:N_fft/2)), 'r', 'LineWidth', 1.5);
grid on;
title('Magnitude Response of FIR Filters');
xlabel('Frequency (rad/sample)');
ylabel('Magnitude (dB)');
legend('Rectangular', 'Hanning', 'Blackman');
axis([0 pi -100 5]);  % Adjust axis for better visibility

%% Assignment #4.3
% FIR Filter Design using Fourier Series Method

N = 71;                           % Filter length
n = -floor(N/2):floor(N/2);       % Symmetric sample indices (centered at 0)
Wc = 0.3 * pi;                    % Cutoff frequency in radians

% Ideal impulse response (low-pass)
hd = (Wc / pi) * sinc((Wc / pi) * n);  % Ideal sinc-based LPF

% Apply different windows
w_rect     = ones(1, N);              % Rectangular window
w_hamming  = hamming(N)';             % Hamming window
w_blackman = blackman(N)';            % Blackman window

% FIR filter coefficients
h_rect     = hd .* w_rect;
h_hamming  = hd .* w_hamming;
h_blackman = hd .* w_blackman;

% Plot and Compare Frequency Responses
fvtool(h_rect, 1, h_hamming, 1, h_blackman, 1);
legend('Rectangular', 'Hamming', 'Blackman');