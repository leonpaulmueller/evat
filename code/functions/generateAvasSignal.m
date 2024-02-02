function out = generateAvasSignal(velocity, coeffs)
%generateAvasSignal - Generates AVAS signal for given velocity vector and coefficient set.
%
% Syntax: out = generateAvasSignal(velocity, coeffs)
%
% Inputs:
%    velocity - velocity input vector
%    coeffs - coefficient struct, generated from getAvasCoeffs.m
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
        velocity (:,1)
        coeffs struct
    end

    %% Check that all required parameters are defined
    assert(isfield(coeffs, 'synthesisMode'), 'Missing setting: synthesisMode' )
    switch coeffs.synthesisMode
        case 'subtractive'
           requiredSettings = {'fs','synthesisMode', 'velocityVector', ...                  
                'windowSize', 'ofb', 'ir', 'gain'};   
        case 'additive'
           if ~isfield(coeffs, 'includeAmplitudeModulation')
                coeffs.includeAmplitudeModulation = false;
           end
           if coeffs.includeAmplitudeModulation
                requiredSettings = {'fs','synthesisMode', 'velocityVector', 'ofb', 'gain', 'nOscillators', 'oscMag', 'oscFreq', ...
                    'oscActiveRange', 'ampModFreq', 'ampModMag'};          
           else
                requiredSettings = {'fs','synthesisMode', 'velocityVector', 'ofb', 'gain', 'nOscillators', 'oscMag', 'oscFreq', ...
                    'oscActiveRange'};       
           end
        case 'sample'
            requiredSettings = {'fs', 'synthesisMode', 'velocityVector', 'ofb', 'gain', 'sample', 'sampleLength', 'sampleMag', 'sampleSpacing'};
        otherwise
            error('Unknown synthesis mode! Choose subtractive, additive or sample.')
    end
    assert(all(isfield(coeffs, requiredSettings)), sprintf('Missing setting %s !\n', requiredSettings{~isfield(coeffs, requiredSettings)}));


    %% Synthesize AVAS Signal
    switch coeffs.synthesisMode
        % Subtractive Synthesis
        case 'subtractive'
            % Divide input velocity in blocks
            [velocityBlocks, velocityRemain] = buffer(velocity, coeffs.windowSize);
            velocitySimulationBlocks = [mean(velocityBlocks, 1), mean(velocityRemain)];
            nSimulationBlocks = size(velocitySimulationBlocks, 2);

            % Initialize filters
            FIR = cell(1,2);
            FIR{1} = dsp.FIRFilter('NumeratorSource','Input port');
            FIR{2} = dsp.FIRFilter('NumeratorSource','Input port');
            
            % Fade Function to blend blocks
            fadeFcn = hann(2*coeffs.windowSize);
            fadeFcn = fadeFcn(1:coeffs.windowSize);
            
            % Loop through all blocks
            out = zeros(coeffs.windowSize, nSimulationBlocks);
            for blockIdx = 1 : nSimulationBlocks
                % Generate block of white noise
                noise = randn(coeffs.windowSize, 1);
                % Apply IR of previous velocity
                if blockIdx > 1
                    out(:,blockIdx) = FIR{1}(noise, coeffs.ir(:,currentVelocityIdx).'); 
                end
                % Get current velocity
                currentVelocityIdx = knnsearch(coeffs.velocityVector', velocitySimulationBlocks(:, blockIdx));
                % Apply IR for current velocity and fade with previous
                out(:,blockIdx) = flip(fadeFcn) .* out(:,blockIdx) + fadeFcn .* FIR{2}(noise, coeffs.ir(:,currentVelocityIdx).');  
            end

            % Transform output to vector, limit to original length
            out = out(:);
            out = out(1:length(velocity));

        % Additive Synthesis
        case 'additive'
            % create oscillators
            for oscIdx = coeffs.nOscillators : -1 : 1
                % Evaluate magnitude and frequency polynominal for input velocity
                tmpMag = polyval(coeffs.oscMag(oscIdx,:), velocity);
                tmpF = polyval(coeffs.oscFreq(oscIdx,:), velocity);
                
                % Limit to fs/2
                tmpMag(tmpF >= coeffs.fs/2) = 0;
                tmpMag(tmpF < 0) = 0;

                % Limit to active velocity range if that range is not
                % covering the entire velocity range
                if coeffs.oscActiveRange(oscIdx, 1) > 1
                    tmpMag(velocity < coeffs.velocityVector(coeffs.oscActiveRange(oscIdx, 1))) = 0;
                end
                if coeffs.oscActiveRange(oscIdx, 2) < length(coeffs.velocityVector)
                    tmpMag(velocity > coeffs.velocityVector(coeffs.oscActiveRange(oscIdx, 2))) = 0;
                end

                % Fade in and out to avoid magnitude jumps
                fadeLength = 1;
                fadeInEndIdx = find(diff(tmpMag == 0) == -1)+1;
                fadeOutStartIdx = find(diff(tmpMag == 0) == 1);

                fadeInVec = zeros(length(tmpMag),1);
                for n = 1 : length(fadeInEndIdx)
                    fadeInStartIdx = max(1, fadeInEndIdx(n) - round(fadeLength * coeffs.fs));
                    fadeInVec(fadeInStartIdx+1 : fadeInEndIdx(n)) =  max(fadeInVec(fadeInStartIdx+1 : fadeInEndIdx(n)), linspace(0, tmpMag(fadeInEndIdx(n)), fadeInEndIdx(n) - fadeInStartIdx )');
                end

                fadeOutVec = zeros(length(tmpMag),1);
                for n = 1 : length(fadeOutStartIdx)
                    fadeOutEndIdx = min(length(tmpMag), fadeOutStartIdx(n) + round(fadeLength * coeffs.fs));
                    fadeOutVec(fadeOutStartIdx(n) : fadeOutEndIdx) =  max(fadeOutVec(fadeOutStartIdx(n) : fadeOutEndIdx), linspace(tmpMag(fadeOutStartIdx(n)), 0, fadeOutEndIdx - fadeOutStartIdx(n)+1 )');
                end
                tmpMag = max([tmpMag, fadeInVec, fadeOutVec], [], 2);
                
                % Apply amplitude modulation
                if coeffs.includeAmplitudeModulation
                    fModVec = polyval(coeffs.ampModFreq(oscIdx,:), velocity);
                    ampModVec = polyval(coeffs.ampModMag(oscIdx,:), velocity);
                    % Create modulator
                    modulator = ampModVec .* cos(2.*pi.*cumsum(fModVec)/coeffs.fs);
                    out(:,oscIdx) = (tmpMag + modulator) .* cos(2.*pi.*cumsum(tmpF)/coeffs.fs);
                else
                    out(:,oscIdx) = tmpMag .* cos(2.*pi.*cumsum(tmpF)/coeffs.fs);
                end
            end

            % Sum up oscillators to output signal
            out = sum(out, 2);
        
        % Sample-based synthesis
        case 'sample'
            out = zeros(length(velocity), 1);
            doneFlag = false;
            currentStartSample = 1;
            while ~doneFlag
                currentEndSample = currentStartSample + coeffs.sampleLength - 1;
                if currentEndSample <= length(out)
                    out(currentStartSample : currentEndSample) = coeffs.sample * (polyval(coeffs.sampleMag, velocity(currentStartSample)) / rms(coeffs.sample));
                    currentStartSample = currentStartSample + round(polyval(coeffs.sampleSpacing, velocity(currentStartSample)) * coeffs.fs);
                else
                    doneFlag = true;
                end
            end
            
        otherwise
            error('Unknown synthesis mode! Choose either subtractive or additive or sample.')
    end

    % Apply Gain compensation 
    if ~isempty(coeffs.gain)
        % Divide signal into octave bands
        release(coeffs.ofb);
        outSignalBands = coeffs.ofb(out);
        
        % Apply gain coefficients
        outSignalBands = outSignalBands .* coeffs.gain;
   
        % remove group delay of octave filter
        groupDelay = round(getGroupDelays(coeffs.ofb));
        outSignalPadded = [outSignalBands;zeros(max(groupDelay),size(outSignalBands, 2))];
        outSignalBandsDelayComp = zeros(size(outSignalBands));
        for bandIdx = 1:size(outSignalBands, 2)
            outSignalBandsDelayComp(:,bandIdx) = outSignalPadded(groupDelay(bandIdx)+1:size(outSignalBands, 1)+groupDelay(bandIdx),bandIdx);
        end

        % Combine bands to output signal
        out = sum(outSignalBandsDelayComp, 2);
    end

end

%------------- END OF CODE --------------