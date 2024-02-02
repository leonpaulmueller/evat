function out = makeBinaural(in, hrtf, coords, blockSize)
%makeBinaural - Apply HRTF to input signal depending on source position to
%obtain binaural output.
%
% Syntax: out = makeBinaural(in, hrtf, coords, blockSize)
%
% Inputs:
%    in - Input signal vector
%    hrtf - HRTF struct
%    coords - Source coordinates from receiver point of view
%    blockSize - Block size for blockwise processing
%    
% Outputs:
%    out - Binaural output signal
%
% Author: Leon Müller
% Email: leon.mueller@chalmers.se
% Website: www.ta.chalmers.se
% January 2024; Last revision: 22/01/2024

%------------- BEGIN CODE --------------
    arguments
        in (:,1)
        hrtf struct
        coords (:,3)
        blockSize (1,1)
    end
    
    % Get spherical source coordinates from POV of receiver
    coords = [coords; repmat(coords(end,:), size(in,1) - size(coords, 1), 1)];
    [az, el] = cart2sph(coords(:,1), coords(:,2), coords(:,3));

    % Receiver is facing street, i.e. azimuth is rotated by 90 degrees
    az = rad2deg(az) - 90;
    az = wrapTo360(az);
    el = rad2deg(el);

    % Divide source position angles in blocks
    az = buffer(az, blockSize);
    el = buffer(el, blockSize);

    % Get mean source position for each block
    az = mean(az, 1)';
    el = mean(el, 1)';
    assert(all(hrtf.SourcePosition(1,:) >= 0))

    % Find nearest neighbor for source angle in HRTF set 
    idx = knnsearch(hrtf.SourcePosition(1:2,:)', [az, el]);

    % Prepare FIR filter objects to apply HRTFs
    FIR = cell(1,2);
    FIR{1} = dsp.FIRFilter('NumeratorSource','Input port');
    FIR{2} = dsp.FIRFilter('NumeratorSource','Input port');
    FIR{3} = dsp.FIRFilter('NumeratorSource','Input port');
    FIR{4} = dsp.FIRFilter('NumeratorSource','Input port');

    % Divide source signal in blocks, allocate space for output
    in = buffer(in, blockSize);
    out = zeros(blockSize, length(az), 2);

    % Prepare fade between old and new HRTF to avoid clicking artifacts
    window = hanning(blockSize*2);
    fadeIn = window(1:blockSize);
    fadeOut = window(blockSize+1:end);

    % Loop over all blocks
    for n = 1 : size(az, 1)
            % Filter current input block with current HRTF
            tmpL = FIR{1}(in(:,n), hrtf.Data_IR(:,1,idx(n)).'); % Left
            tmpR = FIR{2}(in(:,n), hrtf.Data_IR(:,2,idx(n)).'); % Right
        if n > 1
            % Filter current input block with previous HRTF
            tmpLold = FIR{3}(in(:,n), hrtf.Data_IR(:,1,idx(n-1)).');
            tmpRold = FIR{4}(in(:,n), hrtf.Data_IR(:,2,idx(n-1)).');
        else
            % For the first block there is no previous HRTF
            tmpLold = tmpL;
            tmpRold = tmpR;
        end
        % Create output block by fading between previous and current HRTF block
        out(:,n,1) = fadeIn .* tmpL + fadeOut .* tmpLold;
        out(:,n,2) = fadeIn .* tmpR + fadeOut .* tmpRold;
    end
    out = reshape(out, [], 2);
end

%------------- END OF CODE --------------