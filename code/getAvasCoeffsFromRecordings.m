%getAvasCoeffsFromRecordings - Demonstration on how to analyze AVAS or Tire/Road Noise recordings
%to obtain the coefficients required for synthesis with generateAvasSignal.m
%
% Required toolboxes:
%   - Audio Toolbox
%   - DSP System Toolbox
%   - Signal Processing Toolbox
%   - Statistics and Machine Learning Toolbox
%   - Image Processing Toolbox
%
% Tested for Matlab versions >= R2021b
%
% Make sure that the current working directory is the code folder for the
% relative paths to work!
%
% Author: Leon Müller
% Email: leon.mueller@chalmers.se
% Website: www.ta.chalmers.se
% January 2024; Last revision: 01/02/2024

%------------- BEGIN CODE --------------
clear; close all; clc
% Add data and function directory to path
addpath(genpath(['..' filesep 'data']))
addpath(['.' filesep 'functions'])
%% Settings

% Path to recording to be analyzed - If you want to use your own
% recordings, check the provided data to understand the data structure.
recPath = 'vehicleA_AVAS_forwards.mat';

% Decide if settings are defined manually or loaded from mat file
loadSettings = true;
loadSettingsPath = 'vehicleA_forwards_settings.mat';

if loadSettings
    % Load pre-defined analysis settings
    load(loadSettingsPath)
else
    % Define analysis settings manually

    % general settings
    settings.fs = 6e3;                          % Analysis Sampling Rate
    settings.synthesisMode = 'subtractive';     % 'subtractive', 'additive' or 'sample' (see Paper for detailed description)
    settings.velocityResolution = 0.1;          % Resolution of velocity analysis in km/h
    settings.gainCalVelocityRange = [2 25];     % Velocity range that is used for gain compensation
    
    % Substractive synthesis settings, also required for additive synthesis
    if strcmpi(settings.synthesisMode, 'additive') || strcmpi(settings.synthesisMode, 'subtractive')
        settings.fftSize = 2^10;                % FFT size for analysis
        settings.windowSize = 2^8;              % Window size for analyis
        settings.windowOverlap = 0.98;          % Overlap between blocks
        settings.windowType = 'rect';           % 'rect' or 'hanning'
        settings.doMotionFilter = true;         % Smooth magnitude spectrum over velocity
        settings.motionFilterLength = 14;       % Filter length for velocity smoothing (higher number = more smoothing)
        settings.velStdThreshold = 0.05;        % Velocity standard deviation threshold for blocks to be excluded from analysis
    end
    
    % Additive Synthesis Settings
    if strcmpi(settings.synthesisMode, 'additive')
        
        % Peak detection and polynomial fitting
        settings.nOscillators = 27;                 % Total numeber of oscillators (see paper for detailed explanation)
        settings.peakThreshold = -30;               % Threshold for peak detection, relative to maximum
        settings.peakProminence = 2;                % Minimum prominence for peak detection
        settings.peakDistance = 20;                 % Minimum distance between two peaks in Hz
        settings.frequencyPolynomialOrder = 1;      % Polynominal order for Oscillator Frequency
        settings.magnitudePolynomialOrder = 3;      % Polynominal order for Oscillator Magnitude
        settings.ransacSampleSize = 3;              % RANSAC sample size (see help ransac)
        settings.ransacMaxDistance = 25^2;          % RANSAC maximum distance (see help ransac)
        settings.ransacMaxNumTrials = 10000;        % RANSAC maximum number of trials (see help ransac)
        settings.ransacMaxSamplingAttempts = 1000;  % RANSAC maximum number of sampling attempts (see help ransac)
        settings.ransacConfidence = 95;             % RANSAC confidence (see help ransac)
        
        % Amplitude modulation settings
        settings.includeAmplitudeModulation = false; 
        if settings.includeAmplitudeModulation
            settings.modWindowSize = 2^12;          % Window size for modulation analysis
            settings.modWindowOverlap = 0.95;       % Window overlap for modulation analysis
            settings.modFftSize = 2^14;             % FFT size for modulation analysis
            settings.modFmax = 15;                  % Maximum modulation rate
            settings.modThreshold = 0.1;            % Modulation depth detection threshold
            settings.modDetectionBandwith = 0.3;    % Modulation detection frequency bandwith factor
            settings.modPolyOrder = 1;              % Modulation over velocity polynominal order
        end
    end
    
    % Sample based synthesis settings
    if strcmpi(settings.synthesisMode, 'sample')
        settings.referenceRegion = [2.24, 3.23];    % Manually define one period of the repeated sample, in seconds
        settings.corrThreshold = 0.5;               % Exclusion threshold for cross-correlation between reference samples and repetions
        settings.fadeLength = 0.001;                % In/Out fade length for each repetition
    end
end

%% Analyze Signal
[coeffs, referencePressure, generatedPressure] = getAvasCoeffs(recPath, settings);

%% Play recorded and re-synthesized signals
disp('Playing recorded signal')
soundsc(referencePressure, coeffs.fs);
pause(20)
clear sound
pause(1)
disp('Playing generated signal')
soundsc(generatedPressure, coeffs.fs);
pause(20)
clear sound

%% Save coefficients
saveData = false;
if saveData
    save("vehicleX_forwards_coeffs.mat", "coeffs")
end

%------------- END OF CODE --------------