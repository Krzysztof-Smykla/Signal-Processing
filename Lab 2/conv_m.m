function [y, ny] = conv_m(x, nx, h, nh)
% Convolution routine for signal processing
% ------------------------------------------
% [y, ny] = conv_m(x, nx, h, nh)
% [y, ny] = convolution result and its time index
% [x, nx] = first signal and its time indices
% [h, nh] = second signal and its time indices

% Determine the starting and ending indices of the output
nyb = nx(1) + nh(1);                  % Start of output time index
nye = nx(end) + nh(end);             % End of output time index
ny = nyb:nye;                         % Output time index vector

% Perform the convolution
y = conv(x, h);
end

