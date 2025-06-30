% Exercize 1 

% P1a

% generate the signal delta[n] and plot it
clc,clearvars

n = -10: 10;  %values of time domain

delta_n = [zeros(1,10), 1, zeros(1,10)];

%b
figure(1)
stem(n, circshift(delta_n, 2), LineWidth=1); % delta impulse shifted circularly by 2

axis([-10, 10, 0, 1.5]);

title('Unit Sample Function');

xlabel('Time index n');
ylabel('\delta[n]');
%%
%c
figure(2)
u_step = [zeros(1,10), 1, ones(1,10)];
stem(n,u_step,LineWidth=1); % unit step function plot
axis([-10, 10, 0, 1.5]);

title('Unit Step Function');

xlabel('Time index n');
ylabel('U[n]');

%d
figure(3)
stem(n, circshift(-u_step, 3), LineWidth=1); % unit step shifted by v[3, -1]
axis([-10, 10, 0, 1.5]);

title('Unit Step Function 2');

xlabel('Time index n');
ylabel('U[-n - 3]');
%%

% P2a

% generate and plot a discrete-time cosine signal
clc, clearvars

n = 0:40;  % values of the time variable 
w = 0.1*2*pi;  % frequency of the sinusiod.
phi = 0; % phase offset.
A = 1.5; % amplitude
xn = A * cos(w*n - phi); % signal formula

figure(1)
grid("on") % using the grid function
k = stem(n, xn, LineWidth=1);  % defining the sinusoid plot in discrete [n] domain
axis([0, 40, -2, 2]);
grid;
title('Discrete Time Sinusoid');
xlabel('Time index n');
ylabel('x[n]');

%b
l = length(xn) % signal length

fprintf('The length of the signal is %d \n', l )

%c
% METHOD 1

% % Find peaks to determine periodicity
% [pks, locs] = findpeaks(xn, 'MinPeakDistance', 1);
% % using findpeaks() form the Signal Processing toolbox
% 
% 
% % Compute the period (distance between peaks)
% period_estimates = diff(locs);
% 
% % Fundamental period is the most common difference
% N = mode(period_estimates);
% 
% % Display the result
% fprintf('The fundamental period of the discrete-time sinusoid is N = %d samples.\n', N);

%c
% METHOD 2

% The fundamental period can be calculated by finding the difference between the local maxima of the signal. 
% Later count the number of samples from crest to crest

maxYValue = max(xn); 
indices = find(xn == maxYValue); % Get all indices where max value occurs
if length(indices) > 1
    xValues = sum(abs(indices(1) - indices(2))); % Difference of first two occurrences of the maximum
else
    xValues = indices; % If only one occurrence, return the index itself
end
fprintf("The fundamental period is %d samples \n",xValues)

%c
% METHOD 3

% Alterantively, simply count the number of samples present in the signal
% from crest to crest.

%d

% The purpose of the grid() function 
% is to display a grid inside a given
% figure in order to improve visualization. 

%%
clc, clearvars

% P3a

% Use Matlab o generate and plot the discrete-time-signal x[n] = sin(wn)
% for the following values of w:

% -29pi/8,  -3pi/8, -pi/8, pi/8, 3pi/8, 5pi/8, 7pi/8, 9pi/8, 13pi/8,
% 15pi/8, 33pi/8, and 21pi/8 .

n = 0:63;  % discrete-time domain
k_values = [-29, -3, -1, 1, 3, 5, 7, 9, 13, 15, 33, 21];
numPlots = length(k_values);
plotsPerFigure = 4;  % We want a 4x1 grid in each figure

for i = 1:numPlots
    % Open a new figure every time we start a new group of 4 plots
    if mod(i-1, plotsPerFigure) == 0
        figure;
    end

    % Determine subplot index within the current figure (1 to 4)
    subplotIndex = mod(i-1, plotsPerFigure) + 1;
    subplot(plotsPerFigure, 1, subplotIndex);

    % Compute the angular frequency and the corresponding sinusoid
    w = k_values(i) * pi/8;
    x = sin(w * n);

    % Plot the sinusoid using stem
    stem(n, x, 'LineWidth', 1);
    title(sprintf('%d\\pi/8', k_values(i)));
    xlabel('Time index n');
    ylabel('x[n]');
    axis([min(n), max(n), -2, 2]);
    grid on;
end

% After analyzing the visualizations, I concluded that the graph of the
% sinusoid changes its shape every 2𝜋/8 of a rotation.

% The signal repeats (i.e. is periodic) when the argument of the sine
% increases by a multiple of 2𝜋.

% The graphs repeat every 2*k*𝜋 rotations where k = 8.

%%

%P4a

% Consider the Matlab code below which generates a continuous-time complex
% exponential signal and then graphs the real and imaginary parts in one figure and the
% magnitude and phase in another figure.

% generate and plot a continuous-time complex sinusoid

clc, clearvars

t = -4:0.01:4; % values of the time variable.
w = 2.2; % frequency of the sinusoid.

xt = exp(1i*w*t);
xtR = real(xt);
xtI = imag(xt);

figure(1);  % make Fig 1 active
plot(t, xtR, '-b'); % '-b' means 'solid blue line'

 axis([-4, 4, -1, 2]);
 grid on;
 hold on;   % add more curves to the same graph
 plot(t, xtI, '-r');  % solid red line
 title('Real and Imaginary parts');
 xlabel('Time t');
 ylabel('x(t)');
 legend('Re[x(t)]', 'Im[x(t)]');
 hold off;

 mag = abs(xt);
 phase = angle(xt);

 figure(2); % make Fig 2 active
 plot(t, mag, '-g'); % solid green line
 grid on;
 hold on;
 plot(t, phase, '-r'); % solid red line
 title('Magnitute and Phase');
 legend('|x(t)|', 'arg[x(t)]');
 xlabel('Time t');
 ylabel('x(t)');
 hold off;

%%
% % P4a cont.
clc, clearvars

% New signal xt_2

t_2 = 0:0.01:4; % new time domain
w2 = 8; % new fequency
xt_2 = 3 .* exp(- t_2 ./2 ) .* exp(1i*w2*t_2)  % new complex signal

xt_2R = real(xt_2);
xt_2I = imag(xt_2);

figure(3); % Open new figure

plot(t_2, xt_2R, '-b'); % '-b' means 'solid blue line'

 axis([-4, 4, -1, 2]);
 grid on;
 hold on;   % add more curves to the same graph
 plot(t_2, xt_2I, '-r');  % solid red line
 title('Real and Imaginary parts');
 xlabel('Time t_2');
 ylabel('x(t_2)');
 legend('Re[x(t_2)]', 'Im[x(t_2)]');
 hold off;



 mag2 = abs(xt_2);
 phase2 = angle(xt_2);


figure(4);
plot(t_2, mag2, '-g'); % solid green line
 grid on;
 hold on;
 plot(t_2, phase2, '-r'); % solid red line
 title('Magnitute and Phase');
 legend('|x(t_2)|', 'arg[x(t_2)]');
 xlabel('Time t_2');
 ylabel('x(t_2)');
 hold off;

 