function out = applyAirAttenuation(in, coords, fs, humidity, temperature, ambientPressure)
%applyAirAttenuation - Applies air attenuation according to ISO 9613-1:1993
%
% Syntax: out = applyAirAttenuation(in, coords, fs, humidity, temperature, ambientPressure)
%
% Inputs:
%    in - input signal (samples, channels)
%    coords - source coordinates relative to receiver position (samples, xyz, channels)
%    fs - sampling frequency
%    humidity - relative humidity in % (optional, default = 50)
%    temperature - temperature in C (optional, default = 20)
%    ambientPressure - ambient pressure in Pa (optional, default = 101325)
%
% Outputs:
%    out - attenuated signal
%
% Author: Leon Müller
% Email: leon.mueller@chalmers.se
% Website: www.ta.chalmers.se
% January 2024; Last revision: 10/01/2024

%------------- BEGIN CODE --------------

arguments
    in   
    coords (:,3,:)
    fs (1,1)
    humidity (1,1) = 50;
    temperature (1,1) = 20;
    ambientPressure (1,1) = 101325;
end
    % Divide input signal in 1/3 octave bands
    ofb = octaveFilterBank("1/3 octave",fs,'FilterOrder', 4, 'FrequencyRange', [50 fs/2]);
    fc = getCenterFrequencies(ofb);
    in = (ofb(in));

    % Get distance
    d = vecnorm(coords,2,2);  
    
    % Calculate air attenuation
    pa0 = 101325;             % Reference ambient pressure
    T = 273.15 + temperature; % Temperature in Kelvin
    T0 = 293.15;              % Reference temperature
    T01 = 273.16; 
    C = -6.8346 * (T01/T)^1.261 + 4.6151;
    
    % Saturation vapour pressure over water 
    psat = pa0*10^C;
    % Molar concentration of water vapour in %
    h = humidity * psat/ambientPressure;
    % Relaxation frequency of oxygen
    frO = (ambientPressure/pa0) * (24 + 4.04*10^4 * h * (0.02+h)/(0.391+h)); 
    % Relaxation frequency of nitrogen 
    frN = (ambientPressure/pa0)/sqrt(T/T0)*(9+280*h*exp(-4.17*((T/T0)^(-1/3)-1))); 
    
    % Attenuation in dB/m
    att = 8.686*fc.^2.*((1.84*10^(-11)*(pa0/ambientPressure)*sqrt(T/T0))+(T/T0)^(-5/2)*(0.01275*(exp(-2239.1/T)).*... 
        (frO+fc.^2/frO).^(-1)+0.1068*exp(-3352/T).*(frN+fc.^2/frN).^(-1)));
    
    % Multiply attenuation with distance and get attenuation coefficient
    att = 10.^(-(d .* att) / 20);

    % Apply air attenuation to 1/3-octave bands
    out = in .* att;
    
    % Compensate filter bank delay
    groupDelay = round(getGroupDelays(ofb));
    outPadded = [out; zeros(max(groupDelay), size(out, 2), size(out, 3))];
    outDelayComp = zeros(size(out));
    for bandIdx = 1:size(out, 2)
        outDelayComp(:,bandIdx,:) = outPadded(groupDelay(bandIdx)+1:size(out, 1)+groupDelay(bandIdx),bandIdx, :);
    end
    
    % Sum 1/3-octave bands to output signal
    out = squeeze(sum(outDelayComp, 2));
    
end

%------------- END OF CODE --------------