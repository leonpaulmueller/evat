function out = generateOutSignal(srcSignal, irs)
%generateOutSignal - Generates output signal at receiver position by applying moving Green's
%functions to source signal.
%
% Syntax: out = generateOutSignal(srcSignal, irs)
%
% Inputs:
%    srcSignal - Source signal (source samples, channels)
%    irs - Impulse response matrix. Must contain one ir for each source signal sample (ir samples, source samples, channels)
%
% Outputs:
%    out - generated signal
%
% Author: Leon Müller
% Email: leon.mueller@chalmers.se
% Website: www.ta.chalmers.se
% January 2024; Last revision: 11/01/2024

%------------- BEGIN CODE --------------
    arguments
        srcSignal 
        irs 
    end

    disp('--- Generating output signal by applying moving Greens functions ---')
    tt1 = tic;

    % Check that dimensions match
    assert(size(srcSignal,2) == size(irs, 3))
    nSourceSamples = size(srcSignal, 1);
    nIrSamples = size(irs, 1);
    nChannels = size(irs, 3);
    
    % Allocate space for out signal
    out = zeros(nSourceSamples + nIrSamples + 1, 1, nChannels);

    % Pre-calculate irs scaled with source signal
    irs = permute(srcSignal, [3 1 2]) .* irs;

    % Loop over all samples of source signal and add up irs
    for n = 1 : nSourceSamples
        shiftedIdx = n : n + nIrSamples - 1;
        out(shiftedIdx, :, :) = out(shiftedIdx, :, :) + irs(:,n,:);
    end
    
    if nChannels > 1
        out = squeeze(out);
    end

    t = toc(tt1);
    fprintf('--- generateOutSignal duration: %.2f seconds ---\n', t)
end
%------------- END OF CODE --------------