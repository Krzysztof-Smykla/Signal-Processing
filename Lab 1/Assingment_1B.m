% Lab assignment B1
% Difference equation
% y[n] = 1.97y[n-1] - y[n-2]
yn = 1; yn1 = 0;

while i <= y
    i = i+1
    yn2 = yn1; % 0
    yn1 = yn; % 1
    yn = 1.97*yn1-yn2;
end

%%
% A1.1
% y[n] = 1.97*y[n-1]-y[n-2]
clc, clearvars


y = zeros(1, 100);
y(1) = 0
y(2) = 1

for n = 3:100
    y(n) = 1.97*y(n-1)-y(n-2);
end

figure(2)
stem(1:100, y);
xlabel('n')     % Label for x-axis
ylabel('y[n]')  % Label for y-axis
title('Sine Wave Generator using Difference Equation')
axis([0,55,-6, 6])


%%
% A1.2
% y[n] = 1.90*y[n-1]-y[n-2]
clc, clearvars

y = [0, 1, zeros(1, 98)];


for n = 3:100
    y(n) = 1.90*y(n-1) - y(n-2);
end

figure(2)
stem(1:100, y);
xlabel('n')     % Label for x-axis
ylabel('y[n]')  % Label for y-axis
title('Sine Wave Generator using Difference Equation')
axis([0,55,-6, 6])



%%
% A2
% Second order system:
% y[n] = ax[n] + by[n-1] + cy[n-2]


clc; clearvars;

% Define simulation parameters
N = 100;        % Number of time steps (must be equal to length of y)
y = [0,1, zeros(1, 98)]; 

% Modified input parameters a, b, c
cases = [1, -1.5, 0.8; 
         1,  -1.2, 0.8;
         1.45, -1.5, 0.8;
         1,  -1.5, 1.2];

x = zeros(1, 100); % Excitation input x[n]

for i = 1:size(cases,1)
    figure(i); % Create a separate figure for each case
    
    % Extract parameters
    a = cases(i, 1);  % First column
    b = cases(i, 2);  % Second column
    c = cases(i, 3);  % Third column
    
    % Compute system response iteratively
    for n = 3:N  % Start from n=3 since y(1) and y(2) are given
        y(n) = a*x(n) + b*y(n-1) + c*y(n-2);
    end

    % Plot results
    stem(1:N, y,'filled'); % discrete-time plot
    xlabel('n');     
    ylabel('y[n]');  
    title(sprintf('Second Order System: a=%.2f, b=%.2f, c=%.2f', a, b, c));
    axis([0, 33, -6, 6]); % Adjusted axis limits
    grid on;
end



