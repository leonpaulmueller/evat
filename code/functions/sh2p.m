function p = sh2p(X, f,  coords, r0, c, showProgress, memoryLimit)
%sh2p - Get pressure from SH coefficients at evaluation coordinates.
%
% Syntax: p = sh2p(X, f,  coords, r0, c, showProgress)
%
% !! Requires getSH.m !! (https://github.com/polarch/Spherical-Harmonic-Transform/blob/master/getSH.m)
%
% Inputs:
%    X - SH coefficients (f, SH coeffs)
%    f - Frequency vector
%    coords - Cartesian coordinates of pressure evaluation points relative to source position (position nr, xyz, channels)
%    r0 - Radius of evaluation sphere used to determine SH coefficients
%    c - Speed of sound in m/s (default = 343)
%    showProgress - Flag to output progress messages (default = true)
%    memoryLimit - RAM Limit in GB to decide whether to loop over frequencys or to apply one single matrix multiplication (default = 6)
%
% Outputs:
%    p - Pressure at evaluation points
%
% References:
%    [1] Earl G. Williams, Fourier Acoustics, 1999, isbn: 9780127539607
%    [2] Archontis Politis, Microphone array processing for parametric spatial audio techniques, 2016
%        Doctoral Dissertation, Department of Signal Processing and Acoustics, Aalto University, Finland
%
% Author: Leon Müller
% Email: leon.mueller@chalmers.se
% Website: www.ta.chalmers.se
% January 2024; Last revision: 11/01/2024

%------------- BEGIN CODE --------------
    arguments
        X
        f(:,1)
        coords
        r0 (1,1)
        c (1, 1) = 343;
        showProgress = true;
        memoryLimit = 6;
    end

    % Check if getSH.m is on path
    assert(exist('getSH', 'file') == 2, ['getSH.m not found! Please download it from ' ...
        '<a href="https://github.com/polarch/Spherical-Harmonic-Transform/blob/master/getSH.m">' ...
        'https://github.com/polarch/Spherical-Harmonic-Transform/blob/master/getSH.m</a> and add it to the Matlab path!'])

    tt1 = tic;

    % Handle multidimensional input arrays
    reshapedInput = false;
    if ~ismatrix(coords)
        % Reshape evaluation coordinates to 2D matrix
        reshapedInput = true;
        inSize = size(coords);
        coords = permute(coords, [1 3:ndims(coords) 2]);
        coords = reshape(coords, [], 3);      
    end

    % Get attributes of input data
    order = sqrt(size(X, 2)) - 1;
    nPoints = size(coords, 1);
    nFreqs = size(X, 1);
    assert(nFreqs == length(f));

    % Get coordinate angles, assumes source at 0,0
    [az, el, r] = cart2sph(coords(:,1), coords(:,2),coords(:,3));
    
    % Get SH basis functions [2]
    Y = getSH(order, [az, pi/2-el], 'real');

    % The distance of the evaluation points to the center, r, cant be smaller than the reference radius r0
    assert(all(round(r, 3) >= r0), 'r can not be smaller than r0!')

    % Check if any r differs from r0, otherwise we don't need to scale 
    if any(round(r, 3) ~= r0)
        % if r ~= r0 scale SH with hankel functions (See [1] Eq. 6.94)
        p = zeros(nFreqs, nPoints);
        Y = permute(Y, [2 3 1]);

        % Precompute as much as possible
        hnX1 = f.*2.*pi./c.*permute(r, [3 2 1]);
        hnX2 = f.*2.*pi./c.*r0;
        constantSphFactor = sqrt(pi./(2*hnX1)) ./ sqrt(pi./(2*hnX2));   

        % Loop through SH orders
        for n = 0:order
            % Get SH indices for current order
            idx = n^2 + 1 : n^2 + 2*n + 1;  
           
            % Multiply SH coeffs with SH shape function, apply 2nd kind Hankel fraction for distance scaling, sum up over all orders.
            p = p + squeeze(constantSphFactor .* (besselh(n+0.5, 2, hnX1) ./ besselh(n+0.5, 2, hnX2)) .* pagemtimes(X(:, idx), Y(idx, :, :)));
            
            % Show progress
            if showProgress
                switch n
                    case 0 
                        fprintf('sh2p order %i/%i ', n, order)
                    case order
                        fprintf('%i\n', n)
                    otherwise
                        if rem(n, 30) == 0
                            fprintf('\nsh2p order %i/%i ', n, order)
                        else
                            fprintf('%i ', n)
                        end
                end
            end
        end
    else
        % If r = r0 just do matrix multiplication without Hankel scaling
        % Single operation for small matrices
        if numel(X)*size(Y,1)*8/1e9 < memoryLimit
            p = squeeze(sum(X .* permute(Y, [3 2 1]), 2));

        % Loop over frequencies for bigger matrices
        else
            p = zeros(nFreqs, nPoints);
            for fIdx = 1 : nFreqs
                p(fIdx, :) = Y * X(fIdx,:).';

                % Show progress
                if showProgress
                    switch fIdx
                        case 1 
                            fprintf('sh2p f %i/%i ', fIdx, nFreqs)
                        case nFreqs
                            fprintf('%i\n', fIdx)
                        otherwise
                            if rem(fIdx, 30) == 0
                                fprintf('\nsh2p f %i/%i ', fIdx, nFreqs)
                            else
                                fprintf('%i ', fIdx)
                            end
                    end
                end
            end

        end
    end
    
    % Reshape pressure to match original input dimensions
    if reshapedInput
        p = reshape(p, [size(p,1), inSize(1), inSize(3:end)]);   
    end

    t = toc(tt1);
    if showProgress
        fprintf('sh2p duration: %.2f seconds \n', t)
    end

end
%------------- END OF CODE --------------