function X = truncateSH(X, order, doWindow, windowLength)
%truncateSH - Truncates SH coefficients to a lower order.
%
% Syntax: X = truncateSH(X, order, doWindow, windowLength)
%
% Inputs:
%    X - SH coefficients (f, SH coeffs)
%    order - Desired SH order
%    doWindow - Flag to apply hanning window, fading out the highest SH orders (default = false)
%    windowLength - Percentage of SH coefficients to fade out (default = 0.5, i.e. 50%)
%
% Outputs:
%    X - Truncated SH coefficients
%
% Author: Leon Müller
% Email: leon.mueller@chalmers.se
% Website: www.ta.chalmers.se
% January 2024; Last revision: 11/01/2024

%------------- BEGIN CODE --------------
arguments
    X
    order (1,1)
    doWindow = false
    windowLength = 0.5
end
    % Get energy of original SH coefficients
    energy = sum(abs(X).^2, 2);

    % Truncate SH coefficients to lower order 
    oldOrder = sqrt(size(X, 2)) - 1;
    assert(order <= oldOrder, 'New SH order must be <= old SH order!')
    X = X(:,1:(order+1)^2);

    % Apply window 
    if doWindow
        % Generate appropriate Hanning window
        window = hann(round((order+1) * windowLength * 2));
        window = window(ceil((end-1)/2)+1:(end-1));
        window = [ones(order+1 - length(window), 1); window];
        
        % Loop through orders and apply windowing coefficients
        for n = 0:order
            idx = n^2 + 1 : n^2 + 2*n + 1;  % get SH indices for current order
            X(:, idx) = window(n+1) .* X(:, idx); % multiply with window value
        end
            
        % compensate for energy loss due to windowing
        X = X .* energy./sum(abs(X).^2, 2);
    end
end

%------------- END OF CODE --------------