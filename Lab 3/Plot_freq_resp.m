% Lab 3, Task 3

% Define the frequency range
w = linspace(0, pi, 501);           % 501 points from 0 to pi

% Compute the frequency response
H = exp(1j*w) ./ (exp(1j*w) - 0.9);  % Frequency response

% Plot magnitude and phase
figure;

% Magnitude Response
subplot(2,1,1);
plot(w, abs(H), 'b', 'LineWidth', 1.5);
title('Magnitude Response |H(e^{j\omega})|');
xlabel('\omega (rad/sample)');
ylabel('Magnitude');
grid on;

% Phase Response
subplot(2,1,2);
plot(w, angle(H), 'r', 'LineWidth', 1.5);
title('Phase Response ∠H(e^{j\omega})');
xlabel('\omega (rad/sample)');
ylabel('Phase (radians)');
grid on;