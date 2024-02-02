%bem2SH - Demonstration on how to convert the pressure obtained from BEM
%model to spherical harmonic coefficients
%
% External dependencies: 
%   - 'getSH.m' from https://github.com/polarch/Spherical-Harmonic-Transform/blob/master/getSH.m
%
% Required toolboxes:
%   - Statistics and Machine Learning Toolbox
%
% Tested for Matlab versions >= R2021b
%
% Make sure that the current working directory is the code folder for the
% relative paths to work!
%
% Author: Leon Müller
% Email: leon.mueller@chalmers.se
% Website: www.ta.chalmers.se
% January 2024; Last revision: 02/02/2024

%------------- BEGIN CODE --------------
clear; close all; clc
% Add data and function directory to path
addpath(genpath(['..' filesep 'data']))
addpath(['.' filesep 'functions'])

% Settings
order = 64;     % Spherical harmonic order
cutIrs = true;  % Truncate IRs
irLength = 150; % IR truncation length

dataPath = 'vehicleA_3m_5810.mat';

%% Load and pre-process data 
disp('Pre-processing data')

% Check if getSH.m is on path
assert(exist('getSH', 'file') == 2, ['getSH.m not found! Please download it from ' ...
    '<a href="https://github.com/polarch/Spherical-Harmonic-Transform/blob/master/getSH.m">' ...
    'https://github.com/polarch/Spherical-Harmonic-Transform/blob/master/getSH.m</a> and add it to the Matlab path!'])

% Load BEM pressure (frequency domain) on evaluation sphere
tmpData = load(dataPath);
measCoords = tmpData.coords;
p = tmpData.p;
f = tmpData.f;
fs = 2*max(f);
clear tmpData

% Divide pressure by jw to obtain flat frequency spectrum if source would be a
% monopole in free-field
p = p ./ (1i * 2 * pi * f);

% Add DC Bin if not already included in pressure set
if f(1) > 0.1 
    f = [0; f];
    p = [zeros(1, size(p,2)); p];
end

% Get IRs from BEM pressure for plotting and truncating
ir = ifft(ss2ds(p));

% Truncate IR and transform back to single sided pressure
if cutIrs
    % Make irLength even, truncate IR and transform back to frequency domain
    irLength = irLength+rem(irLength,2);
    ir = ir(1:irLength, :);
    p = fft(ir);
    p = p(1:floor(irLength/2)+1, :);
    f =(0:floor(irLength/2))*(fs/irLength).';
end

% Plot Raw BEM pressure and IR
plotIRs(p, f, ir)

%% Fit Spherical Harmonics
disp('Fitting SH coefficients')

% Setup coordinate System based on evaluation points of BEM pressure
[az, el, r0] = cart2sph(measCoords(:,1), measCoords(:,2),measCoords(:,3));
assert(range(r0) < 1e-5, 'Sphere is not centered!')
r0 = round(mean(r0), 3); % assume equal radius to all points

% Get SH basis functions
Y = getSH(order, [az, pi/2-el], 'real');

% Calculate spherical harmonic coefficients by least mean square solution
% regularization or backslash makes no difference here, backslash is faster
% than pinv though
X = (Y \ p.').' ;

%% Plot Error
doErrorPlots = true;
if doErrorPlots
    disp('Calculating error')
    % get pressure at original positions back from SH and check error
    % compared to BEM pressure
    pSH = sh2p(X,f,measCoords,r0);
    plotErrors(p, pSH, f);
end

%% Save Results
saveData = false;
if saveData
    save('shCoeffs.mat', "X", "order", "f", "fs", "r0")
end

%% Functions
function plotIRs(p, f, ir)
    figure('Name', 'Pre-processed BEM pressure')
    
    % Plot magnitude spectrum
    nexttile
    plot(f,20*log10(abs(p(:,1:100:end))));
    xlabel('Frequency in Hz')
    ylabel('Magnitude in dB')
    
    % Plot impulse responses
    nexttile
    plot(ir(:,1:100:end))
    xlabel('n Samples')
    ylabel('Ir')
    
    % Plot maximum magnitude of all impulse responses
    nexttile
    plot(20*log10(max(abs(ir),[],2)./max(abs(ir), [], "all" )))
    xlabel('n samples')
    ylabel('max(Ir) dB')

    % Plot mean magnitude of all pressure points
    nexttile
    plot(f,20*log10(mean(abs(p), 2)));
    xlabel('Frequency in Hz')
    ylabel('Mean Magnitude')
end


function [rmsError, MAC] = plotErrors(p, pEst, f)
    % Calculate normalized RMS error and modal assurance criterion
    MAC = abs(diag(p*pEst')).^2 ./ (diag(p*p') .* diag(pEst*pEst'));
    rmsError = 20*log10(rms(abs(p-pEst), 2)./rms(p, 2));

    figure('Name', 'SH error compared to BEM pressure')
    tiledlayout flow
    
    % Plot RMS error
    nexttile
    semilogx(f, rmsError.')
    xlabel('Frequency in Hz')
    ylabel('RMS Error in dB')
    xlim([50 3000])
    
    % Plot modal assurance criterion
    nexttile
    semilogx(f, MAC.')
    ylim([0 1])
    xlabel('Frequency in Hz')
    ylabel('MAC')
    xlim([50 3000])
end

%------------- END OF CODE --------------