function [coeffs, referencePressure, generatedPressure, referenceVelocity] = getAvasCoeffs(fileName, coeffs)
%getAvasCoeffs - Get AVAS synthesis coefficients from reference recording.
%
% Syntax: [coeffs, referencePressure, generatedPressure, referenceVelocity] = getAvasCoeffs(fileName, coeffs)
%
% Inputs:
%    fileName - path to AVAS recording. See provided avas measurements for required file structure 
%    coeffs - Analysis settings. See getAvasCoeffsFromRecordings.m for detailed description of required settings.
%
% Outputs:
%    coeffs - Determined AVAS synthesis coefficients
%    referencePressure - Recorded pressure signal
%    generatedPressure - Re-synthesized pressure signal
%    referenceVelocity - Velocity vector for recorded AVAS signal
%
% Author: Leon Müller
% Email: leon.mueller@chalmers.se
% Website: www.ta.chalmers.se
% January 2024; Last revision: 01/02/2024

%------------- BEGIN CODE --------------
    arguments
        fileName char
        coeffs struct
    end
    
    %% Check that all relevant settings are defined
    coeffs.srcFileName = fileName;
    assert(isfield(coeffs, 'synthesisMode'), 'Missing setting: synthesisMode' )
    switch coeffs.synthesisMode
        case 'subtractive'
           requiredSettings = {'fs','synthesisMode' ,'velocityResolution','gainCalVelocityRange','velStdThreshold','fftSize', ...                  
        'windowSize' ,'windowOverlap','windowType','doMotionFilter','motionFilterLength'};   
        case 'additive'
           if ~isfield(coeffs, 'includeAmplitudeModulation')
                coeffs.includeAmplitudeModulation = false;
           end
           if coeffs.includeAmplitudeModulation
                requiredSettings = {'fs','synthesisMode' ,'velocityResolution','gainCalVelocityRange','velStdThreshold','fftSize', ...                  
                'windowSize' ,'windowOverlap','windowType','doMotionFilter','motionFilterLength','nOscillators', ...              
                'peakThreshold','peakProminence','peakDistance','frequencyPolynomialOrder','magnitudePolynomialOrder', ...
                'ransacSampleSize','ransacMaxDistance','ransacMaxNumTrials','ransacMaxSamplingAttempts','ransacConfidence', ...          
                'includeAmplitudeModulation','modWindowSize','modWindowOverlap','modFftSize','modFmax','modThreshold', ...              
                'modDetectionBandwith','modPolyOrder'};          
           else
                requiredSettings = {'fs','synthesisMode' ,'velocityResolution','gainCalVelocityRange','velStdThreshold','fftSize', ...                  
                'windowSize' ,'windowOverlap','windowType','doMotionFilter','motionFilterLength','nOscillators', ...              
                'peakThreshold','peakProminence','peakDistance','frequencyPolynomialOrder','magnitudePolynomialOrder', ...
                'ransacSampleSize','ransacMaxDistance','ransacMaxNumTrials','ransacMaxSamplingAttempts','ransacConfidence'};     
           end
        case 'sample'
            requiredSettings = {'fs', 'synthesisMode', 'velocityResolution', 'gainCalVelocityRange', 'referenceRegion','corrThreshold','fadeLength'};
        otherwise
            error('Unknown synthesis mode! Choose subtractive, additive or sample.')
    end
    assert(all(isfield(coeffs, requiredSettings)), sprintf('Missing setting %s !\n', requiredSettings{~isfield(coeffs, requiredSettings)}));
    
    if ~isfield(coeffs, 'doOnlyFRF')
        coeffs.doOnlyFRF = false;
    end

    %% Load reference data
    tmpFile = load(fileName);
    fsRecording = tmpFile.fs;
    pressureCell = tmpFile.p;
    vehicleVelocityCell = tmpFile.v;
    
    % Combine all individual measurements to one long recording
    referencePressure = cat(1,tmpFile.p{:});
    referenceVelocity = cat(1,tmpFile.v{:});

    % Resample to desired fs
    referencePressure = resample(referencePressure, coeffs.fs, fsRecording); 
    referenceVelocity = resample(referenceVelocity, coeffs.fs, fsRecording); 
    
    % Determine maximum measured velocity
    if ~isfield(coeffs, 'maxVelocity')
        coeffs.maxVelocity = max(referenceVelocity);
    end
    
    %% Analyze reference signals
    switch coeffs.synthesisMode

        % Estimate FRF for additive and subtractive synthesis modes
        case {'additive', 'subtractive'}
            % Divide recording in blocks for each velocity step
            [pressureBlocksSortedFRF, coeffs.velocityVector] = divideInBlocks(pressureCell, vehicleVelocityCell, fsRecording,  coeffs.windowSize, ...
                round(coeffs.windowSize * coeffs.windowOverlap), coeffs);
           
            % Apply window to blocks
            switch coeffs.windowType
                case 'hanning'
                    window = hanning(coeffs.windowSize);
                    window = window ./ sqrt( sum(window .^ 2) / coeffs.windowSize );
                case 'rect'
                    window = 1;
                otherwise 
                    error('Unknown window type! Choose either rect or hanning')
            end
    
            % Calculate averaged magnitude spectra FRF
            frf = zeros(coeffs.fftSize, length(pressureBlocksSortedFRF));
            for cellIdx = 1 : length(pressureBlocksSortedFRF)
                if isempty(pressureBlocksSortedFRF{cellIdx})
                    frf(:,cellIdx) = nan;
                else
                    frf(:,cellIdx) = mean( abs( fft( pressureBlocksSortedFRF{cellIdx} .* window, coeffs.fftSize, 1 ) ./ coeffs.fftSize ), 2 );
                end
            end
            clear pressureBlocksSortedFRF
    
            % Remove velocities above highest estimated FRF
            fullIdx = find(~isnan(frf(1,:)));
            frf = frf(:,1:fullIdx(end));
            coeffs.velocityVector = coeffs.velocityVector(1:fullIdx(end));
            coeffs.nVelocityValues = length(coeffs.velocityVector);
            coeffs.maxVelocity = coeffs.velocityVector(end);
            
            % If empty set FRF for 0 km/h to zero
            if isnan(frf(1,1))
                frf(:,1) = zeros(coeffs.fftSize, 1);
                fullIdx = [1, fullIdx];
            end
            
            % Interpolate FRF for empty entries
            frfMagInterp = interp1(coeffs.velocityVector(fullIdx),abs(frf(:,fullIdx)).',coeffs.velocityVector, "linear").';
            
            % Apply motion filter in order to avoid jumps between velocity steps
            motionFilter = fspecial('motion',coeffs.motionFilterLength,0) ;
            if coeffs.doMotionFilter
                coeffs.frf = imfilter(abs(frfMagInterp), motionFilter) ;
            else
                coeffs.frf = frfMagInterp;
            end
            
            % Get minimum phase impulse response and frequency vector
            coeffs.ir = getMinPhaseFilt(ifft(coeffs.frf));
            coeffs.fFRF = linspace(0, coeffs.fs, coeffs.fftSize);

            if coeffs.doOnlyFRF
                return
            end
    
            % Plot resulting FRF
            figure('Name', 'Estimated FRF')  
            srf = surf(coeffs.velocityVector, coeffs.fFRF, 20*log10(abs(coeffs.frf ./ 2e-5)), 'LineStyle','none');
            ylim([10 coeffs.fs/2])
            xlim([min(coeffs.velocityVector) max(coeffs.velocityVector)])
            xlabel('Velocity in km/h')
            ylabel('Frequenccy in Hz')
            cb = colorbar;
            ylabel(cb, 'Level in dB SPL')
            view([0 90])
            
            % Analyze FRF for additive synthesis
            if strcmpi(coeffs.synthesisMode, 'additive')
                % Update plot for additive synthesis
                ax = gca;
                set(srf, 'FaceAlpha', 0.9)
                view([-15 80])
                ax.ZAxis.Visible = 'off';
                ax.ZGrid = 'off';
                ax.Color = 'none';
                grid off
                hold on
        
                % Allocate loop variables
                peakFrequencies = cell(coeffs.nVelocityValues, 1);
                peakMagnitudes = cell(coeffs.nVelocityValues, 1);
                tmpVelVec = cell(coeffs.nVelocityValues, 1);
                
                % Loop through all velocity steps and find peaks
                for velocityIdx = 1:coeffs.nVelocityValues
                    % Find maximum of current velocity step
                    maxval = max(abs(coeffs.frf(1:coeffs.fftSize/2,velocityIdx)));
        
                    % Transform magnitudes to log relative to current maximum and find peaks
                    [peakMagnitudes{velocityIdx}, peakFrequencies{velocityIdx}] = findpeaks(20*log10(abs(coeffs.frf(1:coeffs.fftSize/2,velocityIdx)./maxval)), ...
                        coeffs.fFRF(1:coeffs.fftSize/2), 'MinPeakProminence', coeffs.peakProminence, 'MinPeakHeight',coeffs.peakThreshold, 'MinPeakDistance', coeffs.peakDistance);
        
                    % Transform magnitudes back to linear
                    peakMagnitudes{velocityIdx} = maxval * 10.^(peakMagnitudes{velocityIdx}./20);
                    tmpVelVec{velocityIdx} = repmat(coeffs.velocityVector(velocityIdx), size(peakMagnitudes{velocityIdx}, 1), 1);
                end
                    
                % Combine cells of different velocity steps to one matrix
                peakMagnitudes = cat(1, peakMagnitudes{:});
                tmpVelVec = cat(1, tmpVelVec{:});
                peakFrequencies = cat(2, peakFrequencies{:})';

                % Mark peaks in plot
                scatter3(tmpVelVec, peakFrequencies,20*log10(abs(peakMagnitudes./ 2e-5)), 'xg')
                
                % Combine frequencies and magnitudes with corresponding velocity 
                peakFrequencies = [tmpVelVec, peakFrequencies];
                peakMagnitudes = [tmpVelVec, peakMagnitudes];
        
                % Define fit and evaluation functions for RANSAC
                fitLineFcn = @(xy) polyfit(xy(:,1),xy(:,2),coeffs.frequencyPolynomialOrder ); 
                evalLineFcn = @(P, xy) sum((xy(:, 2) - polyval(P, xy(:,1))).^2,2); 
        
                % Allocate loop variables
                coeffs.oscFreq = zeros(coeffs.nOscillators, coeffs.frequencyPolynomialOrder+1);
                coeffs.oscMag = zeros(coeffs.nOscillators, coeffs.magnitudePolynomialOrder+1);
                coeffs.oscActiveRange = zeros(coeffs.nOscillators, 2);
      
                % Loop through all oscillators and perform RANSAC for frequency and normal polyfit for magnitude
                for n = 1 : coeffs.nOscillators
                    
                    % Do RANSAC fitting for frequency
                    [coeffs.oscFreq(n,:), inlierIdx] = ransac(peakFrequencies,fitLineFcn,evalLineFcn, ...
                      coeffs.ransacSampleSize,coeffs.ransacMaxDistance,"MaxNumTrials", coeffs.ransacMaxNumTrials, ...
                      "MaxSamplingAttempts",coeffs.ransacMaxSamplingAttempts, "Confidence",coeffs.ransacConfidence);
                    
                    % Determine velocity range in which current oscillator is active
                    coeffs.oscActiveRange(n,:) = [find(ismember(coeffs.velocityVector, peakFrequencies(inlierIdx, 1)), 1, 'first'), ...
                                                    find(ismember(coeffs.velocityVector, peakFrequencies(inlierIdx, 1)), 1, 'last')];
            
                    % Use maximum of all inliers for each velocity step to define magnitude
                    currentMag = peakMagnitudes(inlierIdx,:);
                    currentMag = sortrows(currentMag);
                    [~,maxIdx] = unique(currentMag(:,1),'last');
                    currentMag = currentMag(maxIdx, :);
        
                    % Interpolate missing magnitude values within active velocity range
                    currentMag = interp1(currentMag(:,1), currentMag(:,2), coeffs.velocityVector(coeffs.oscActiveRange(n,1):coeffs.oscActiveRange(n,2)));
                    
                    % Fit polynominal for magnitude of current oscillatore
                    coeffs.oscMag(n,:) = polyfit(coeffs.velocityVector(coeffs.oscActiveRange(n,1):coeffs.oscActiveRange(n,2)),currentMag,coeffs.magnitudePolynomialOrder);
                    
                    % Remove points that have been used by current RANSAC iteration from dataset
                    peakFrequencies(inlierIdx, :) = [];
                    peakMagnitudes(inlierIdx, :) = [];
        
                    % Plot polynominal for current oscillator
                    plotVel = coeffs.velocityVector(coeffs.oscActiveRange(n,1) :coeffs.oscActiveRange(n,2));
                    plotFreq = polyval(coeffs.oscFreq(n,:),plotVel);
                    plotMag = 20*log10(abs(polyval(coeffs.oscMag(n,:),plotVel)./2e-5));
                    plot3(plotVel,plotFreq,plotMag, 'r', 'LineWidth', 3)
        
                end
                hold off
    
                % Analyze Amplitude modulation
                if coeffs.includeAmplitudeModulation
    
                    % Divide recordings in blocks
                    [pressureBlocksSortedMod, coeffs.velocityVector] = divideInBlocks(pressureCell, vehicleVelocityCell, ...
                        fsRecording, coeffs.modWindowSize, round(coeffs.modWindowSize * coeffs.modWindowOverlap), coeffs);
    
                    % Create modulation frequency vector
                    fVecMod = linspace(0, coeffs.fs, coeffs.modFftSize);
                    [~, fMaxModIdx] = min(abs(fVecMod - coeffs.modFmax));
                    fVecMod = fVecMod(1:fMaxModIdx);
                
                    % Prepare lowpass filter for envelope detection
                    lpFiltEnv = firpm(20,[0 0.03 0.1 1],[1 1 0 0]);
                
                    % Get averaged envelope for each velocity step
                    for velIdx = length(coeffs.velocityVector) : -1 : 1    
                        if ~isempty(pressureBlocksSortedMod{velIdx})
    
                            % Bandpass filter pressure signal for each oscillator frequency
                            clear bandFilteredBlock
                            for oscIdx = coeffs.nOscillators : -1 : 1
                                tmpfc = polyval(coeffs.oscFreq(oscIdx,:), coeffs.velocityVector(velIdx));
                                bandFilteredBlock(:, oscIdx, :) = permute(bandpass(pressureBlocksSortedMod{velIdx}, ...
                                    [tmpfc * (1 - coeffs.modDetectionBandwith), tmpfc * (1 + coeffs.modDetectionBandwith)], ...
                                    coeffs.fs, ImpulseResponse="fir"), [1 3 2]);
                            end
                            
                            % Get envelope of BP filtered signals 
                            envelopeBlock = abs(hilbert(bandFilteredBlock));
    
                            % Lowpass filter envelope and remove DC offset
                            for n = 1 : size(envelopeBlock, 3)
                                envelopeBlock(:,:,n) = filter(lpFiltEnv, 1, envelopeBlock(:,:,n));
                                envelopeBlock(:,:,n) = envelopeBlock(:,:,n) - rms(envelopeBlock(:,:,n)); % remove DC offset
                            end
                            
                            % Average envelope magnitude in frequency domain
                            tmpModSpec =  mean(abs(fft(envelopeBlock,coeffs.modFftSize)/coeffs.modFftSize), 3); 
                            modSpec(:, velIdx, :) = tmpModSpec(1:fMaxModIdx, :,:);
                        else
                            modSpec(:, velIdx, :) = nan(fMaxModIdx, 1, coeffs.nOscillators);                                                                
                        end
                    end
                    clear pressureBlocksSortedMod
    
                    % Find empty velocities
                    emptyIdx = isnan(modSpec(1,:,1));
                    fullIdx = find(~emptyIdx);
                
                    % Remove empty velocities at end
                    modSpec = modSpec(:,1:fullIdx(end), :);
                
                    % If empty set frf for 0 km/h to zero
                    if fullIdx(1) > 1
                        modSpec(:,1,:) = zeros(size(modSpec, 1), 1, size(modSpec, 3));
                    end
                
                    % Truncate velocity vector to valid range
                    coeffs.velocityVector = coeffs.velocityVector(1:fullIdx(end));
                    coeffs.maxVelocity = coeffs.velocityVector(end);
                    
                    % Interpolate empty velocities
                    modSpec = permute(interp1(coeffs.velocityVector(fullIdx),permute(modSpec(:,fullIdx,:), [2 1 3]),coeffs.velocityVector, "linear"), [2 1 3]); 
                    
                    % Ignore DC bin
                    modSpec(1,:) = 0;
                
                    % Plot modulation spectrum for each oscillator
                    figure('Name', 'Modulation Spectrum')
                    tiledlayout flow
                    for n = 1 : size(modSpec, 3)
                        nexttile
                        surf(coeffs.velocityVector, fVecMod, abs(modSpec(:,:,n)), 'LineStyle','none') %  20*log10(abs(modSpec(:,:,n))+0.1)
                        cb = colorbar;
                        ylabel(cb, 'Modulation Amplitude')
                        view([0 90])
                        title(['Modulation Spectrum Oscillator ' num2str(n)])
                        xlabel('Velocity in km/h')
                        ylabel('Modulation Frequency')
                    end
                
                    % Find frequency bin with strongest modulation for each oscillator and velocity
                    [ampMod, fModIdx] = max(modSpec);
                    fMod = squeeze(fVecMod(fModIdx)); 
                    ampMod = squeeze(ampMod);
                    
                    % Set modulation amplitudes below threshold to zero
                    ampMod(ampMod < coeffs.modThreshold) = 0;
                    fMod(ampMod < coeffs.modThreshold) = 0;
                
                    % Fit polynominals to modulation frequency and amplitude
                    for modIdx = size(modSpec, 3) : -1 : 1
                        coeffs.ampModMag(modIdx,:) = polyfit(coeffs.velocityVector,ampMod(:,modIdx),coeffs.modPolyOrder);
                        coeffs.ampModFreq(modIdx,:) = polyfit(coeffs.velocityVector,fMod(:,modIdx),coeffs.modPolyOrder);
                    end
                
                    % Update velocity range of previous results
                    coeffs.nVelocityValues = min(size(coeffs.frf, 2), length(coeffs.velocityVector) );
                    coeffs.frf = coeffs.frf(:,1:coeffs.nVelocityValues);
                    coeffs.ir = coeffs.ir(:,1:coeffs.nVelocityValues);
                    coeffs.oscActiveRange(:,2) = min(coeffs.oscActiveRange(:,2), length(coeffs.velocityVector));
                end                       
            end

        % Sample-based synthesis
        case 'sample'
    
            figure('Name', 'Recorded AVAS Signal'); 
            t = linspace(0, length(referencePressure)/coeffs.fs, length(referencePressure));
            plot(t,referencePressure);
            hold on    
            
            % Cut out reference signal according to user input
            refSignal = referencePressure(round(coeffs.referenceRegion(1) * coeffs.fs) : round(coeffs.referenceRegion(2) * coeffs.fs));
            coeffs.sampleLength = round(diff(coeffs.referenceRegion*coeffs.fs));
            coeffs.velocityVector = 0:coeffs.velocityResolution:coeffs.maxVelocity; % Range of velocity values
            
            % Analyze correlation between reference signal and full recording
            [c,lags] = xcorr(referencePressure, refSignal); 
            c = c(find(lags == 0) : end);
            [~, maxCorrIdx] = findpeaks(c, "MinPeakDistance", coeffs.sampleLength);
            
            % Ignore first and last block to avoid inclomplete samples
            maxCorrIdx = maxCorrIdx(2 : end-1);
            
            % Plot detected sample regions
            xline(t(maxCorrIdx), 'r');
            xline(t(maxCorrIdx) + t(coeffs.sampleLength), 'g');
            
            % Cut out signal segments that have highest correlation with reference        
            nBlocks = length(maxCorrIdx);
            pBlocks = zeros(coeffs.sampleLength, nBlocks);
            vBlocks = zeros(1, nBlocks);
            for n = 1 : nBlocks
                pBlocks(:, n) = referencePressure(maxCorrIdx(n) : maxCorrIdx(n) + coeffs.sampleLength - 1);
                vBlocks(n) = mean(referenceVelocity(maxCorrIdx(n) : maxCorrIdx(n) + coeffs.sampleLength - 1));
            end
            
            % Calculate distance between samples vs vehicle velocity
            blockDistance = diff(maxCorrIdx)/coeffs.fs;
            
            % Fit polynominal to block distance and magnitude
            coeffs.sampleSpacing =  polyfit(vBlocks(1:end-1), blockDistance, 2);
            coeffs.sampleMag = polyfit(vBlocks, rms(pBlocks), 2);
            
            % Plot results of polyfit
            figure('Name', 'Sample Spacing and Magnitude Polynomials')
            tiledlayout('flow')
            nexttile
            scatter(vBlocks(1:end-1), blockDistance)
            hold on
            plot(sort(unique(vBlocks)), polyval(coeffs.sampleSpacing, sort(unique(vBlocks))), 'r')
            xlabel('Velocity in km/h')
            ylabel('Sample Distance in s')
            nexttile
            scatter(vBlocks, rms(pBlocks))
            hold on
            plot(sort(unique(vBlocks)), polyval(coeffs.sampleMag, sort(unique(vBlocks))), 'r')
            xlabel('Velocity in km/h')
            ylabel('RMS Magnitude in Pa')
            
            % Analyze correlation between blocks
            r = corr(pBlocks);
            r(logical(eye(nBlocks))) = nan;
    
            % Calculate mean correlation for each block
            r = tanh(mean(atanh(r), 'omitnan'));
            
            % Sort out blocks with low correlation to rest of blocks
            pBlocks(:, r < coeffs.corrThreshold) = [];
            r(r < coeffs.corrThreshold) = [];
    
            % Calculate overall mean correlation
            rTotal = tanh(mean(atanh(r), 'omitnan'));
            nBlocks = size(pBlocks, 2);
            fprintf('Averaged %i blocks with a total correlation of %.3f \n', nBlocks, rTotal)
            
            % Calculate time domain average
            coeffs.sample = mean(pBlocks, 2);
            
            % Fade in and out
            coeffs.fadeLength = round(coeffs.fadeLength * coeffs.fs);
            fadeWindow = hanning(coeffs.fadeLength * 2);
            coeffs.sample(1 : coeffs.fadeLength) = coeffs.sample(1 : coeffs.fadeLength) .* fadeWindow(1:coeffs.fadeLength);
            coeffs.sample(end-coeffs.fadeLength+1 : end) = coeffs.sample(end-coeffs.fadeLength+1 : end) .* fadeWindow(coeffs.fadeLength+1 : end);
    
            % Plot averaged signal
            figure('Name', 'Averaged Signal')
            tiledlayout flow
            nexttile
            plot(coeffs.sample)
            nexttile
            semilogx(20*log10(abs(fft(coeffs.sample))))
            nexttile
            N = 256;
            spectrogram(coeffs.sample,hanning(N),N-1,N,coeffs.fs, "Yaxis")
    end
    
    
    %% Re-Synthesize reference signal
    coeffs.gain = [];
    coeffs.ofb = [];
    generatedPressure = generateAvasSignal(referenceVelocity, coeffs);
    
    %% Compensate overall gain offset between reference and re-synthesized signal    
    coeffs.ofb = octaveFilterBank("1/3 octave",coeffs.fs,'FilterOrder', 6, 'FrequencyRange', [22 coeffs.fs/2]);
    coeffs.fc = getCenterFrequencies(coeffs.ofb);
    
    release(coeffs.ofb)
    outSignalBands = coeffs.ofb(generatedPressure);
    pressureBands = coeffs.ofb(referencePressure);
    gainCompIdx = (referenceVelocity > max(coeffs.gainCalVelocityRange(1), min(coeffs.velocityVector))) & ...
        (referenceVelocity < min(coeffs.gainCalVelocityRange(2), max(coeffs.velocityVector)));
    coeffs.gain = rms(pressureBands(gainCompIdx,:)) ./ rms(outSignalBands(gainCompIdx, :));
    
    % Re-generate reference signal with gain calibration
    generatedPressure = generateAvasSignal(referenceVelocity, coeffs);
    
    
    %% Compare reference and generated signal
    figure('Name', 'Recorded and Generated Signal')
    windowSizePlot = 2^11;
    tiledlayout flow
    nexttile
    spectrogram(referencePressure/2e-5, hanning(windowSizePlot), 0, windowSizePlot, coeffs.fs,'yaxis')
    ylim([0.01 3]);
    % Older matlab versions don't know clim
    if exist('clim', 'file') == 2
        clim([0 100]);
    elseif exist('caxis', 'file') == 2
        caxis([0 100]);
    end
    title('Recorded Signal')
    nexttile
    spectrogram(generatedPressure/2e-5, hanning(windowSizePlot), 0, windowSizePlot, coeffs.fs,'yaxis')
    title('Generated Signal')
    ylim([0.01 3]);
    if exist('clim', 'file') == 2
        clim([0 100]);
    elseif exist('caxis', 'file') == 2
        caxis([0 100]);
    end
    nexttile
    pwelch(referencePressure/2e-5,windowSizePlot,[],[],coeffs.fs)
    hold on
    pwelch(generatedPressure/2e-5,windowSizePlot,[],[],coeffs.fs)
    legend('Recorded Signal', 'Generated Signal')
    
    levelDiff = 20*log10(rms(highpass(referencePressure, 50, coeffs.fs))/rms(highpass(generatedPressure, 50, coeffs.fs)));
    fprintf('getAvasCoeffs - Overall level difference between recording and synthesis: %.2f dB\n', levelDiff)

    % Arrange figures so that all are visible
    tmpFig = findall(groot,'type','figure');
    if ~isempty(tmpFig)
        for n = 1 : length(tmpFig)
            tmpFig(n).Position = tmpFig(n).Position + [(n-1)*100, 0, 0, 0];
        end
    end
end

%% Divide reference recordings into blocks
function [pressureBlocksSortedFRF, velocityVector] = divideInBlocks(pressureCell, vehicleVelocityCell, fsRecording, windowSize, overlap, coeffs)

        for cellIdx = length(pressureCell) : -1 : 1
            if length(pressureCell{cellIdx}) > 1
                % Resample recordings to desired fs
                tmpPressure = resample(pressureCell{cellIdx},coeffs.fs, fsRecording); 
                tmpVehicleVelocity = resample(vehicleVelocityCell{cellIdx},coeffs.fs, fsRecording);

                % Divide pressure and velocity in blocks, process velocity
                [pressureBlocks, ~] = buffer(tmpPressure,windowSize, overlap, "nodelay");
                [velocityBlocks, ~] = buffer(tmpVehicleVelocity,windowSize, overlap, "nodelay");
                velocityBlocksMean = mean(velocityBlocks, 1);
                velocityBlockIdx = std(velocityBlocks, 1) <= coeffs.velStdThreshold;
                velocityBlocksMean = velocityBlocksMean(:, velocityBlockIdx);
                pressureBlocks = pressureBlocks(:, velocityBlockIdx);
                
                % Get discrete velocity
                velocityBlocksMean = round(velocityBlocksMean ./ coeffs.velocityResolution) * coeffs.velocityResolution;
                velocityVector = 0:coeffs.velocityResolution:coeffs.maxVelocity; % Range of velocity values
                nVelocityValues = length(velocityVector);
                
                % Sort pressure blocks according to velocity
                pressureBlocksSorted = cell(nVelocityValues, 1);
                for velIdx = 1 : nVelocityValues
                    pressureBlocksSorted{velIdx} = pressureBlocks(:, velocityBlocksMean == velocityVector(velIdx));
                end
                pressureBlocksSortedFRFCell{cellIdx} = pressureBlocksSorted;

            end
        end

        % Combine blocks from different recordings to one cell array
        pressureBlocksSortedFRFCell = cat(2, pressureBlocksSortedFRFCell{:});
        for cellIdx = size(pressureBlocksSortedFRFCell, 1) : -1 : 1
            pressureBlocksSortedFRF{cellIdx, 1} = cat(2, pressureBlocksSortedFRFCell{cellIdx, :});
        end
end

%------------- END OF CODE --------------
