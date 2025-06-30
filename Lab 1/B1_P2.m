% Second order system:
% y[n] = ax[n] + by[n-1] + cy[n-2]

% Assignment 1B, p.2

% Use Matlab to generate the discrete feedback systems
% with different parameters a, b and c, plot results.

% The system described in the block diagram is described by the following
% difference eqation
% y[n] = a · x[n] − a · b · y[n − 1] − a · c · y[n − 2]

% Defining parameters
n = 0:60;
y = [0 ones(1, length(n)-1)]; % First value is 0, rest are ones
x = ones(1, length(n));

% V1
a = 1;
b = -1.5;
c = 0.8;
for k = 3:length(n)
y(k) = a*x(k) - a*b*y(k-1) - a*c*y(k-2); % corrected equation based on diagram
end

% Plot
figure(1)
stem(n, y, 'filled');
title(['Difference Equation: a = ', num2str(a), ', b = ', num2str(b), ', c = ', num2str(c),]);
xlabel('Time index n');
ylabel('y[n]');
grid on;

%V2
a = 1;
b = -1.5;
c = 0.6;
for k = 3:length(n)
y(k) = a*x(k) - a*b*y(k-1) - a*c*y(k-2); % corrected equation based on diagram
end

% Plot
figure(2)
stem(n, y, 'filled');
title(['Difference Equation: a = ', num2str(a), ', b = ', num2str(b), ', c = ', num2str(c),]);
xlabel('Time index n');
ylabel('y[n]');
grid on;

%V3
a = 1;
b = -1.5;
c = 0.5;
for k = 3:length(n)
y(k) = a*x(k) - a*b*y(k-1) - a*c*y(k-2); % corrected equation based on diagram
end

% Plot
figure(3)
stem(n, y, 'filled');
title(['Difference Equation: a = ', num2str(a), ', b = ', num2str(b), ', c = ', num2str(c),]);
xlabel('Time index n');
ylabel('y[n]');
grid on;