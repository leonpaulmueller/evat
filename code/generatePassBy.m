%generatePassBy - Generate a simple binaural electric vehicle passage
%
% This script is a basic demonstration of the Electric Vehicle Auralization
% Toolbox (EVAT). It generates a static binaural pass-by with either constant 
% vehicle velocity or a linear acceleration. This code could be extended 
% to, for example, use pre-recorded velocity signals, include image sources or 
% implement a dynamic, head-tracked rendering.
%
% Requires pre-calculated directivities as well as AVAS and tire/road noise
% synthesis coefficients which are provided in data directory or, alternatively, 
% can be generated using getAvasCoeffsFromRecordings.m. 
%
% External dependencies: 
%   - 'getSH.m' from https://github.com/polarch/Spherical-Harmonic-Transform/blob/master/getSH.m
%   - 'HMSII.sofa' from https://sofacoustics.org/data/database/thk/HMSII.sofa
% 
% Required toolboxes:
%   - Audio Toolbox
%   - DSP System Toolbox
%   - Signal Processing Toolbox
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

%% Settings

%---- General Settings ----
c = 343;                    % Speed of sound in m/s
fsOut = 48e3;               % Sample rate of output signal

%---- Pass-by geometry and vehicle velocity ----
velocity = [5, 25, 10];     % Vehicle velocity in km/h. 
                            % Single number for constant velocity
                            % Multiple values for linear acceleration
                            % between values
dur = 20;                   % Duration of pass by
streetDistance = 8;         % Distance between street and receiver
leftToRight = true;         % Movement direction from receiver POV
backwards = false;          % Driving direction from vehicle POV
zPos = 1.8;                 % Height of receiver
tireLevel = -12;            % Level of tire noise in dB relative to AVAS
ambienceLevel = 0;          % Level of ambience noise relative to AVAS and Tire noise

%---- Spatial Resolution ----
shOrderAvas = 32;                   % Spherical Harmonic Order for AVAS
shOrderTires = 8;                   % Spherical Harmonic Order for Tire/Road Noise
spatialDownsamplingFactorAvas = 4;  % Spatial Downsampling Factor for AVAS
spatialDownsamplingFactorTires = 10;% Spatial Downsampling Factor for Tires
% (Higher number reduces spatial resolution and computational demand, a 
% factor of 1 means no spatial downsampling and is the safest choice)

%---- Input data: Change for other vehicle types -----
avasCoeffsPath = 'vehicleA_forwards_coeffs.mat'; % AVAS Synthesis Coefficients
tireCoeffsPath = 'vehicleA_tire_coeffs.mat';     % Tire/Road Noise Synthesis Coefficients
avasSHPath = 'vehicleA_N64Fs6kHz.mat';           % AVAS Spherical Harmonic Directivity
tireSHPath = 'tire_N16Fs16kHz.mat';              % Tire/Road Noise Spherical Harmonic Directivity
ambiencePath = 'vehicleA_ambience.mat';          % Ambience Noise
hrtfPath = 'HMSII.sofa';                         % HRTF Dataset


%% Check that external data and dependencies are available

% Check currrent working directory
if ~strcmpi(flip(extractBefore(flip(pwd), filesep)), 'code')
    warning('Change working directory to code folder, otherwise relative paths might not work!')
end

% Add function directory to path
addpath(['.' filesep 'functions'])

% Add data directory to path
addpath(genpath(['..' filesep 'data']))

% Check if HRTF is available
assert(exist(hrtfPath, 'file') == 2, ['HMSII.sofa HRTF not found! Please download it from ' ...
        '<a href="https://sofacoustics.org/data/database/thk/HMSII.sofa">' ...
        'https://sofacoustics.org/data/database/thk/HMSII.sofa</a> and add it to the Matlab path! Alternatively, change the code to use your own HRTF!'])

% Check if getSH.m is on path
assert(exist('getSH', 'file') == 2, ['getSH.m not found! Please download it from ' ...
        '<a href="https://github.com/polarch/Spherical-Harmonic-Transform/blob/master/getSH.m">' ...
        'https://github.com/polarch/Spherical-Harmonic-Transform/blob/master/getSH.m</a> and add it to the Matlab path!'])

% Start timer
ttTotal = tic;


%% Generate velocity vector and source coordinates at output fs

% Generate velocity vector: Here we could also use recorded velocity data
if length(velocity) == 1
    % Constant velocity
    velVec = repmat(velocity, fsOut * dur, 1);
else
    % Changing velocity with linear acceleration
    velVec = interp1(0 : (dur-1/fsOut)/(length(velocity)-1) : dur-1/fsOut, ...
        velocity, 0 : 1/fsOut : dur-1/fsOut);
end

% Define coordinates of receiver from moving vehicle point of view - this seems
% weird but makes sense in order to get the radiation directivity right!

% Get  x coordinates by integrating velocity vector
coordsSrcPOV(:,1) = -cumsum(velVec / 3.6 / fsOut);

% Make sure that vehicle reaches receiver X position after half of duration
coordsSrcPOV(:,1) = coordsSrcPOV(:,1) - coordsSrcPOV(round(size(coordsSrcPOV, 1)/2),1);
coordsSrcPOV(:,2) = abs(streetDistance);
coordsSrcPOV(:,3) = abs(zPos);

% Define source coordinates from receiver point of view - these coordinates
% are required for the binaural rendering.
coordsRecPOV = coordsSrcPOV .* [1 1 -1];

% Change orientation depending on desired driving direction
if leftToRight
   coordsSrcPOV(:,2) = -coordsSrcPOV(:,2);
   coordsRecPOV(:,1) = -coordsRecPOV(:,1);
end

if backwards
   coordsSrcPOV(:,2) = -coordsSrcPOV(:,2);
   coordsSrcPOV(:,1) = -coordsSrcPOV(:,1);
end


%% Generate AVAS Signal
disp('--- Generating AVAS Signal ---')

% Load AVAS synthesis coefficients and directivity
shAvas = load(avasSHPath);
coeffsAvas = load(avasCoeffsPath, "coeffs");
coeffsAvas = coeffsAvas.coeffs;

% fs is the minimum of source signal and directivity fs
fsAvas = min(shAvas.fs, coeffsAvas.fs);

% Downsample velocity vector to AVAS source signal fs
velVecAvas = resample(velVec, coeffsAvas.fs, fsOut);

% Synthesize AVAS source signal
srcSignalAvas = generateAvasSignal(velVecAvas, coeffsAvas);
srcSignalAvas = resample(srcSignalAvas, fsAvas, coeffsAvas.fs); 

% Downsample source coordinates to AVAS fs
coordsAvas = resample(coordsSrcPOV, fsAvas, fsOut);

% Calculate AVAS Irs describing propagation from source to receiver position
irsAvas = getIRs(shAvas, coordsAvas, fsAvas, shOrderAvas, c, spatialDownsamplingFactorAvas);

% Apply air attenuation
srcSignalAvas = applyAirAttenuation(srcSignalAvas, coordsAvas, fsAvas);

% Generate Avas out signal
outAvas = generateOutSignal(srcSignalAvas, irsAvas);

% Resample to desired output sampling rate
outAvasResampled = resample(outAvas, fsOut, fsAvas);

% Make binaural
sofaInfo = ncinfo(hrtfPath);
for n = 1 : length(sofaInfo.Variables)
    hrtf.(strrep(sofaInfo.Variables(n).Name, '.', '_'))= ncread(hrtfPath, sofaInfo.Variables(n).Name);
end
outAvasResampledBinaural = makeBinaural(outAvasResampled, hrtf, coordsRecPOV, 512);


%% Generate tire noise
disp('--- Generating Tire Noise ---')

% Load tire/road noise synthesis coefficients and directivity
shTires = load(tireSHPath);
coeffsTires = load(tireCoeffsPath, "coeffs");
coeffsTires = coeffsTires.coeffs;

% Tire fs is minimum of synthesis fs and directivity fs
fsTires = min(coeffsTires.fs, shTires.fs);

% Downsample velocity vector to Tire/road source signal fs
velVecTires = resample(velVec, coeffsTires.fs, fsOut);

% Generate tire/road noise for all four tires
for tireIdx = 4 : -1 : 1
    srcSignalTires(:,tireIdx) = generateAvasSignal(velVecTires, coeffsTires);
end

% Resample source signal to tire fs
srcSignalTires = resample(srcSignalTires, fsTires, coeffsTires.fs);

% Downsample source coordinates to Tire fs
tmpCoordsTires = resample(coordsSrcPOV, fsTires, fsOut);

% Offset and rotation matrixes for all four tires relative to center of car
tireOffset = cat(3, [1.5 0.9 0], [-1.5 0.9 0], [1.5 -0.9 0], [-1.5 -0.9 0]);
tireRotation = cat(3, [1 1 1], [1 1 1], [1 -1 1], [1 -1 1]);

% Tire coordinates from source and receiver POV
coordsTiresSrcPOV = tireRotation .* resample(coordsSrcPOV, fsTires, fsOut) + tireOffset;
coordsTiresRecPOV = coordsRecPOV + tireOffset;

% Apply air attenuation to all four tire source signals
srcSignalTires = applyAirAttenuation(srcSignalTires, coordsTiresSrcPOV, fsTires);

% Loop through all tires, get IRs and generate output signals
for tireIdx = 4 : -1 : 1
    fprintf('Tire %i / 4\n', 5-tireIdx)
    % Get IRs describing propagation from source to receiver positions 
    irsTires = getIRs(shTires, coordsTiresSrcPOV(:,:,tireIdx), fsTires, shOrderTires, c, spatialDownsamplingFactorTires);
    % Generate output signal
    outTires{tireIdx} = generateOutSignal(srcSignalTires(:, tireIdx), irsTires);
end

% Combine tire signals to single matrix
maxTireSignalSamples = cellfun(@(x) size(x, 1), outTires, 'UniformOutput', false);
maxTireSignalSamples = max([maxTireSignalSamples{:}]);
outTires = cellfun(@(x) [x; zeros(maxTireSignalSamples - size(x, 1), 1)], outTires, 'UniformOutput', false);
outTires = cat(2, outTires{:});

% Resample to desired output fs 
outTiresResampled = resample(outTires, fsOut, fsTires);
% outTiresResampled = cat(3, outTiresResampled, outTiresResampled);

% Make binaural
for n = 4 : -1 : 1
    outTiresResampledBinaural(:,n,:) = permute(makeBinaural(outTiresResampled(:,n), hrtf, coordsTiresRecPOV(:,:,n), 512), [1 3 2]);
end

% Add up all tires
outTiresResampledBinaural = squeeze(sum(outTiresResampledBinaural, 2));

% Normalize tires to specific level below AVAS
outTiresResampledBinaural = 10^(tireLevel/20) * rms(outAvasResampledBinaural, "all") * outTiresResampledBinaural / rms(outTiresResampledBinaural, "all");

% Zeropadd AVAS and Tire signals to common length
nOutputSamples = max(length(outTiresResampledBinaural), length(outAvasResampledBinaural));
outTiresResampledBinaural = [outTiresResampledBinaural; zeros(nOutputSamples - length(outTiresResampledBinaural), size(outTiresResampledBinaural, 2))];
outAvasResampledBinaural = [outAvasResampledBinaural; zeros(nOutputSamples - length(outAvasResampledBinaural), size(outTiresResampledBinaural, 2))];

% Add up AVAS and tire/road noise
outAvasTires = outTiresResampledBinaural + outAvasResampledBinaural;

%% Add ambience noise
ambience = load(ambiencePath);

% Pick random noise from cell array that is longer than required duration
tmpIdx = find(ambience.dur >= dur+1);
noise = ambience.binaural{tmpIdx(randi(length(tmpIdx)))};

% At noise with desired level relative to AVAS and Tire noise RMS values
out = outAvasTires + 10^(ambienceLevel/20) * rms(outAvasTires, "all") * noise(1:nOutputSamples, :) / rms(noise, "all");

tTotal = toc(ttTotal);
fprintf('--- Done! Total duration: %.2f seconds ---\n', tTotal)

%% Play result
soundsc(out, fsOut)


%% Plot Spectrogram
figure('Name', 'Auralization Result')
N = 2096;
[sSpec,fSpec,tSpec] = spectrogram(out(:,1), hanning(N), round(N*0.8), N*2, fsOut);
surf(tSpec, fSpec, 20*log10(abs(sSpec)/2e-5), 'linestyle', 'none')
view([0 90])
cb = colorbar;
ylabel(cb, 'dB SPL')
set(gca, "YScale", "log")
xlabel('Time in s')
ylabel('Frequency in Hz')
axis tight
ylim([20 5000])

%------------- END OF CODE --------------
