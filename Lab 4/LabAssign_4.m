
% Plotting h(n)
n = -10:20;  % Sample indices for n in range (-20,20)

% Define the ideal impulse response 
h = @(n) 0.3 .* sinc(0.3 .* n);  

ht = zeros(size(n));  % Initialize truncated impulse response

% Defining the truncation logic
for i = 1:length(n)
    if n(i) >= -7 && n(i) <= 7 % if i in range(-7,7)
        ht(i) = h(n(i));  % Call the impulse response with value n(i)
    else
        ht(i) = 0;
    end
end

% Shifted evaluation: define shifted index
n_shifted = n - 7;           % Desired argument: n - 7
ht_shifted = zeros(size(n)); % Create empty matrix

% Lab assign 4.1
% Collect valid shifted samples
vals = [];  % Initialize empty array to store the values

% Safely sample ht using shifted values
for i = 1:length(n)
    idx = find(n == n_shifted(i));  % Find index of (n - 7) in n
    if ~isempty(idx)                % Check if index of ht = ht_shifted
        ht_shifted(i) = ht(idx);
        vals(end+1) = ht(idx);      % Append to list of values
    end
end

% Display the list of truncated values after shift
disp('Truncated impulse response values evaluated at (n - 7):');
disp(vals);

% Plot the truncated impulse response h(n-7)
figure(1);
stem(n, ht_shifted, 'b'); grid on;
title("Truncated impulse response at h(n - 7)");
xlabel('n'); ylabel('h(n - 7)');

% Desinging the FIR filter using fft
% Using the fft methd, design a length -15 FIR lowpass filter to
% approximate an ideal lowpass filter with (omega) W = 0.3pi rad.
N = 15;                 % Length of FIR filter
Wc = 0.3 * pi;          % Cutoff frequency in radians

% Frequency sampling points
k = 0:N-1;
omega = 2 * pi * k / N;

% Desired frequency response: H_d[k] = 1 for |omega| <= 0.3*pi
Hd = double(omega <= Wc | omega >= 2*pi - Wc);  % Symmetric around pi

% Inverse FFT to get time-domain impulse response
h_fir = real(ifft(Hd));  % Impulse response in time domain

% Shift for linear-phase symmetric FIR
h_fir = fftshift(h_fir);  % Center the impulse response

% Display and plot the result
disp('FIR Filter Coefficients:');
disp(h_fir);

% Plot impulse response
figure(2);
stem(0:N-1, h_fir, 'filled');
title('FIR Filter Impulse Response (via FFT method)');
xlabel('n'); ylabel('h[n]');
grid on;

% Optional: frequency response of designed FIR
figure(3);
freqz(h_fir, 1, 512);  % Magnitude and phase response
title('Frequency Response of FIR Filter');

%%
% Assignment 4.2