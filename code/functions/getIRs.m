function irs = getIRs(shCoeffs, coords, fs, order, c, spatialDownSamplingFactor)
%getIRs - Extrapolate SH directivity to obtain Green's functions for all source positions.
%
% Syntax: irs = getIRs(shCoeffs, coords, fs, order, c, spatialDownSamplingFactor)
%
% Inputs:
%    shCoeffs - Spherical harmonic coefficients for sound source directivity (f, sh)
%    coords - Coordinates of sound source relative to receiver position (samples, xyz)
%    fs - Sampling rate
%    order - Desired SH order
%    c - Speed of sound
%    spatialDownSamplingFactor - Factor to downsample spatial resolution of 
%     coordinates to reduce computational demand (>=1, Default = 1)
%    
% Outputs:
%    irs - Green's functions for all source positions
%
% Author: Leon Müller
% Email: leon.mueller@chalmers.se
% Website: www.ta.chalmers.se
% January 2024; Last revision: 11/01/2024

%------------- BEGIN CODE --------------
    arguments
        shCoeffs struct
        coords (:,3)
        fs (1,1)
        order (1,1)
        c (1,1) = 343
        spatialDownSamplingFactor (1,1) = 1;
    end

    tt1 = tic;
    disp('--- starting getIRs ---')

    % Downsample spatial resolution of coordinates to reduce computational demand
    if spatialDownSamplingFactor > 1
        disp('applying spatial downsampling')
        downSamplingVec = 1 : spatialDownSamplingFactor : size(coords, 1);
        highSamplingVec = 1 : size(coords, 1);
        if downSamplingVec(end) ~= size(coords, 1)
            downSamplingVec = [downSamplingVec, size(coords, 1)];
        end
        coords = coords(downSamplingVec, :, :);
    end

    % Truncate SH coefficients to selected order
    if shCoeffs.order ~= order
        shCoeffs.X = truncateSH(shCoeffs.X, order, false);   
    end
    shCoeffs.f = shCoeffs.f(:);

    % Get transfer function for each step
    p = sh2p(shCoeffs.X, shCoeffs.f,  coords, shCoeffs.r0, c);
        
    % Add f=0 to frequency vector if not already included
    if ~ismember(0, shCoeffs.f)
        shCoeffs.f = [0;shCoeffs.f];
        p = [zeros(1, size(p,2), size(p,3)); p];
    else
        % Make sure DC bin is zero
        p(1, :, :) = 0;
    end

    % Get distance between source and receiver to compensate runtime if neccessary
    [~, ~, r] = cart2sph(coords(:,1,:), coords(:,2,:), coords(:,3,:));
    
    disp('compensating runtime delay')
    % Compensate runtime delay:
    % The problem here is that the runtime delay is typically too large to fit into the
    % relatively short IRs which, due to the periodicity of the fourier transform, 
    % leads to a circshift effect. Therefore, we remove the runtime delay in frequency
    % domain, go to time domain, append some zeros and add the original delay
    % in frequency domain. 
    
    % Remove runtime delay from IR
    wDelayBackwards = 2 * pi * shCoeffs.f; 
    delayTermBackwards = exp(-1i .* wDelayBackwards .* -(permute(r - shCoeffs.r0, [2 1 3])/c));
   
    % Multiply pressure with delay compensation term, make double sided, do ifft
    irsTmp = ifft(ss2ds(p .* delayTermBackwards));
    clear delayTermBackwards p

    % Resample to fs
    if fs ~= shCoeffs.fs
        irsTmp = resample(irsTmp, fs, shCoeffs.fs);
    end
    
    % Zeropadd IRs to allow long runtime delay
    maxDelaySamples = ceil(max(fs*(r - shCoeffs.r0)/c, [], 'all'));
    irsTmp = [irsTmp; zeros(maxDelaySamples+1, size(irsTmp, 2), size(irsTmp, 3))]; % was +size(irsTmp, 1)
    
    % Make sure that length is odd so that the spectrum includes no nyquist bin
    if rem(size(irsTmp, 1), 2) == 0
        irsTmp = [irsTmp; zeros(1, size(irsTmp, 2), size(irsTmp, 3))];
    end
    
    % Construct double sided delay factor to re-apply to IRs
    wDelayForwards = 2 * pi * linspace(0, fs/2, ceil(size(irsTmp,1)/2))';
    delayTermForwards = exp(-1i .* wDelayForwards .* (permute(r - shCoeffs.r0, [2 1 3])/c));
    delayTermForwards = ss2ds(delayTermForwards, 1, false); 
    
    % Go back to frequency domain, re-apply original delay, go back to time domain
    irs = ifft(fft(irsTmp) .* delayTermForwards);

    % Interpolate between IRs to compensate for reduced spatial resolution
    if spatialDownSamplingFactor > 1
        irs = permute(interp1(downSamplingVec, permute(irs, [2 1 3]), highSamplingVec), [2 1 3]);
    end
    t = toc(tt1);
    fprintf('--- getIRs duration: %.2f seconds ---\n', t)
   
end

%------------- END OF CODE --------------