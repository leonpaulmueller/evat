function h = getMinPhaseFilt(in)
%getMinPhaseFilt - Get minimum-phase representation of input IR.
%Method according to Matlab 'rceps' documentation
%
% Syntax: h = getMinPhaseFilt(in)
%
% Inputs:
%    in - Input impulse response (samples, channels)
%
% Outputs:
%    h - Minimum-phase representation of input IR
%
% Author: Leon Müller
% Email: leon.mueller@chalmers.se
% Website: www.ta.chalmers.se
% January 2024; Last revision: 11/01/2024

%------------- BEGIN CODE --------------

    blockSize = size(in, 1);
    y = real(ifft(log(abs(fft(in)))));
    w = [1; 2*ones(blockSize/2-1, 1); ones(1-rem(blockSize, 2), 1); zeros(blockSize/2-1, 1)];
    h = real(ifft(exp(fft(w.*y))));
end

%------------- END OF CODE --------------