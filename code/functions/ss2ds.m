function spec = ss2ds(spec, dim, includesNyquist, scale)
%ss2ds - Converts single-sided to double-sided spectrum. Assumes that input
%spectrum includes DC bin.
%
% Syntax: spec = ss2ds(spec, dim, includesNyquist, scale)
%
% Inputs:
%    spec - Input single-sided spectrum
%    dim - Operating dimension (default = 1)
%    includesNyquist - Flag if source spectrum includes Nyquist bin (default = true)
%    scale - Scale so that the energy in double-sided spectrum is the
%    same as in single-sided spectrum (default = false)
%
% Outputs:
%    spec - Double-sided spectrum
%
% Author: Leon Müller
% Email: leon.mueller@chalmers.se
% Website: www.ta.chalmers.se
% January 2024; Last revision: 11/01/2024

%------------- BEGIN CODE --------------
    arguments
        spec
        dim = 1
        includesNyquist = true
        scale = false;
    end
    
    % Shift samples to first dimension
    if dim ~= 1
        spec = shiftdim(spec,dim-1);
    end

    reshapedInput = false;
    if ~ismatrix(spec)
        % Reshape input to 2D matrix
        reshapedInput = true;
        inSize = size(spec);      
        spec = reshape(spec, size(spec,1), []);
    end

    % Construct double sided spectrum
    if includesNyquist
        % DC,f1,..,fn,Ny,-fn,..,-f1
        if scale
            spec = [real(spec(1,:)); spec(2:end-1,:)/2; real(spec(end,:)); conj(flipud(spec(2:end-1,:)/2))];
        else
            spec = [real(spec(1,:)); spec(2:end-1,:); real(spec(end,:)); conj(flipud(spec(2:end-1,:)))];
        end
    else
        % DC,f1,..,fn,-fn,..,-f1
        if scale
            spec = cat(1, spec, conj(flipud(spec(2:end,:))))/2;
            spec(1,:) = 2*real(spec(1,:));
            
        else
            spec = cat(1, spec, conj(flipud(spec(2:end,:))));
            spec(1,:) = real(spec(1,:));
        end
    end
    
    if reshapedInput
        % Reshape spectrum to original input dimensions 
        spec = reshape(spec, [size(spec,1), inSize(2:end)]);
    end

    if dim ~= 1
        spec = shiftdim(spec,-(dim-1));
    end

end

%------------- END OF CODE --------------