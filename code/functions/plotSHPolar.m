function fig = plotSHPolar(X, f, r0, opts)
%%plotSHPolar - Plot polar representation of SH directivity
%
% Syntax: fig = plotSHPolar(X, f, r0, 'Name', 'Value')
%
% Inputs:
%    X - SH Coefficient set
%    f - Frequency vector
%    r0 - Reference radius for SH coefficients
%
% Outputs:
%    fig - Figure handle
%
% Name/Value Properties:
%    'octaveRange' - see bandwith option in octaveFilter documentation, default = "1 Octave".
%    'octaveFilterOrder' - default = 4
%    'fRange' - frequency range for plot, default = [200 max(f)]
%    'minPlotLevel' - plot level threshold for the three subplots, default = [-60, -60, -60]
%    'angularResolution' - angular resoution in degree, default = 1
%    'horizontalPlaneElevationLimits' - Elevation limits for horizontal plane, default = [0 30]
%    'medianPlaneAzimuthLimits' - Azimuth limits for median plane, default = [-45 45]
%    'frontalPlaneAzimuthLimits' - Azimuth limits for frontal plane, default = [30 150]
%    'lineWidth' - plot line width, default = 1.4
%    'imgPath' - path of image for overlay, cell with three separate paths, one for each plane
%    'imgAlpha' - transparency of overlay image, default = 0.8
%    'imgPos' - position of overlay image (3 images x 4 position coords)
%    'rTickRes' - resolution of polar ticks for each plot, default = [10 10 10]
%    'normalization' - "total" for normalization to total maximum, "band" for normalization to maximum of each frequency band, default = "total"
%    'lineStyle' - plot line style for each band, default = {'-', '--', '-', '-.'}
%    'plotHorMedFront' - flag to enable different subplots, default = [1 1 1]
%    'figureName' - Name of figure, default = "SH Directivity"
%    'subPlotTitles' - titles for subplots, default = {'Horizontal Plane', 'Median Plane', 'Frontal Plane'}
%
% Author: Leon Müller
% Email: leon.mueller@chalmers.se
% Website: www.ta.chalmers.se
% February 2024; Last revision: 07/02/2024
%------------- BEGIN CODE --------------
arguments
    X
    f (:,1)
    r0 (1,1) = 1
    opts.octaveRange = "1 Octave"
    opts.fRange = [200,max(f)]
    opts.minPlotLevel = [-60, -60, -60]
    opts.octaveFilterOrder = 4;
    opts.angularResolution = 1;
    opts.horizontalPlaneElevationLimits = [0 30]
    opts.medianPlaneAzimuthLimits = [-45 45]
    opts.frontalPlaneAzimuthLimits = [30 150]
    opts.lineWidth = 1.4
    opts.axFontSize = 14
    opts.legendFontSize = 14
    opts.imgAlpha = 0.8;
    opts.imgPath cell = []
    opts.imgPos (:,4) = [0.37 1.86 0 0; 0.115 0.75 0 0; 0.104 -2.9 15 1]
    opts.rTickRes = [10 10 10];
    opts.normalization {mustBeMember(opts.normalization, ["total", "band"])} = "total"
    opts.lineStyle = {'-', '--', '-', '-.'};
    opts.plotHorMedFront = [1 1 1];
    opts.figureName string = "SH Directivity";
    opts.subPlotTitles = {'Horizontal Plane', 'Median Plane', 'Frontal Plane'};
end

    %% Get pressure from SH coefficients

    % Create equally spaced azimuth and elevation vectors
    azVec = linspace(0, 2*pi, round(360/opts.angularResolution))';
    elVec = linspace(0, pi/2, round(90/opts.angularResolution))';

    % Create evaluation grid
    [azMesh,elMesh,rMesh] = meshgrid(azVec,elVec,r0);
    [xGrid, yGrid, zGrid] = sph2cart(azMesh,elMesh,rMesh);

    % Get pressure on evaluation grid
    fprintf('------ plotSHPolar ------\n')
    p = sh2p(X, f,  [xGrid(:), yGrid(:), zGrid(:)], r0);

    % Add DC bin if not included
    if f(1) > 0.1
        p = [zeros(1, size(p, 2)); p];
    end

    %% Convert pressure to IRs and apply octave filter bank
    irs = ifft(ss2ds(p));
    ofb = octaveFilterBank(opts.octaveRange,2*max(f),'FrequencyRange',opts.fRange, 'FilterOrder', opts.octaveFilterOrder);  
    irsOct = ofb(irs);

    % Get center frequencies and round to even numbers
    fc = getCenterFrequencies(ofb);
    fc(fc <= 500) = roundn(fc(fc <= 500), 1);
    fc(fc > 500) = roundn(fc(fc > 500), 2);
    
    % Get RMS of each filtered IR band and reshape to evaluation grid
    rmsPressure = squeeze(rms(irsOct, 1));
    rmsPressure = reshape(rmsPressure, length(fc), size(azMesh, 1), size(azMesh, 2));

    %% Get level on different evaluation planes
    switch opts.normalization
        case 'total'
            normDim = 'all';
        case 'band'
            normDim = 2;
    end
    
    % Horizontal plane
    if opts.plotHorMedFront(1)
        horizontalPlaneElIdx = elMesh(:,1) >= deg2rad(opts.horizontalPlaneElevationLimits(1)) & ...
            elMesh(:,1) <= deg2rad(opts.horizontalPlaneElevationLimits(2));
        rmsPressureHor = squeeze(sum(rmsPressure(:,horizontalPlaneElIdx,:), 2));
        rmsLevelHorizontal = 20*log10(abs(rmsPressureHor) ./ max(abs(rmsPressureHor), [], normDim));
        rmsLevelHorizontal(rmsLevelHorizontal < opts.minPlotLevel(1)) = opts.minPlotLevel(1);
    end

    % Median plane
    if opts.plotHorMedFront(2) 
        medianPlaneRightAzIdx = azMesh(1,:) >= deg2rad(wrapTo360(opts.medianPlaneAzimuthLimits(1))) | azMesh(1,:) <= deg2rad(wrapTo360(opts.medianPlaneAzimuthLimits(2)));
        medianPlaneLeftAzIdx = azMesh(1,:) >= deg2rad(wrapTo360(opts.medianPlaneAzimuthLimits(1)-180)) & azMesh(1,:) <= deg2rad(wrapTo360(opts.medianPlaneAzimuthLimits(2)-180));
        % remove double value for 0 and 360
        if any(round(rad2deg(azMesh(1,medianPlaneRightAzIdx))) == 0) && any(round(rad2deg(azMesh(1,medianPlaneRightAzIdx))) == 360)
            medianPlaneRightAzIdx( round(rad2deg(azMesh(1,:))) == 360 ) = 0;
        end
        rmsPressureMedianRight = sum(rmsPressure(:, :, medianPlaneRightAzIdx), 3);
        rmsPressureMedianLeft = flip(sum(rmsPressure(:, :, medianPlaneLeftAzIdx), 3), 2); % flip elevation of left half
        % take mean of 90 degree value
        rmsPressureMedian =  [rmsPressureMedianRight(:,1:end-1), mean([rmsPressureMedianRight(:,end), rmsPressureMedianLeft(:,1)], 2), rmsPressureMedianLeft(:,2:end)]; 
        rmsLevelMedian = 20*log10(abs(rmsPressureMedian) ./ max(abs(rmsPressureMedian), [], normDim));
        rmsLevelMedian(rmsLevelMedian < opts.minPlotLevel(2)) = opts.minPlotLevel(2);
    end

    % Frontal plane
    if opts.plotHorMedFront(3)
        frontalPlaneRightAzIdx = azMesh(1,:) >= deg2rad(wrapTo360(opts.frontalPlaneAzimuthLimits(1))) & azMesh(1,:) <= deg2rad(wrapTo360(opts.frontalPlaneAzimuthLimits(2)));
        frontalPlaneLeftAzIdx = azMesh(1,:) >= deg2rad(wrapTo360(opts.frontalPlaneAzimuthLimits(1)+180)) & azMesh(1,:) <= deg2rad(wrapTo360(opts.frontalPlaneAzimuthLimits(2)+180));
        rmsPressureFrontalRight = sum(rmsPressure(:, :, frontalPlaneRightAzIdx), 3);
        rmsPressureFrontalLeft = flip(sum(rmsPressure(:, :, frontalPlaneLeftAzIdx), 3), 2); % flip elevation of left half
        rmsPressureFrontal =  [rmsPressureFrontalRight(:,1:end-1), mean([rmsPressureFrontalRight(:,end), rmsPressureFrontalLeft(:,1)], 2), rmsPressureFrontalLeft(:,2:end)];
        rmsLevelFrontal = 20*log10(abs(rmsPressureFrontal) ./ max(abs(rmsPressureFrontal), [], normDim));
        rmsLevelFrontal(rmsLevelFrontal < opts.minPlotLevel(3)) = opts.minPlotLevel(3);
    end
    
    %% Do plot

    % Create figure
    nPlots = sum(opts.plotHorMedFront);
    switch nPlots
        case 1
            figPos = [177 671 600 376];
        case 2
            figPos = [177 671 1000 376];
        case 3
            figPos = [177 671 1403 376];
    end
    fig = figure('Name', opts.figureName, 'Position', figPos, 'Color',[1 1 1]);

    plotIdx = 1;

    if opts.plotHorMedFront(1)
        % Horizontal Plane
        subplot(1,nPlots,plotIdx)
        % Polar Plot
        polarplot(azMesh(1,:), rmsLevelHorizontal, 'LineWidth', opts.lineWidth)
        rlim([opts.minPlotLevel(1) 0])
        rticks(opts.minPlotLevel(1):opts.rTickRes(1):0);
        rl = rticklabels;
        rl(:,2) = {'dB'};
        rl = join(rl);
        rl(1:2) = {''};
        rticklabels(rl);
        axHor = gca;
        axHor.FontSize = opts.axFontSize;
        axHor.Color = [1 1 1 0];
        axHor.ThetaColor = [0 0 0];
        axHor.RColor = [0 0 0];
        axHor.LineStyleOrder = opts.lineStyle;
        if ~isMATLABReleaseOlderThan("R2023a")
            axHor.LineStyleCyclingMethod = "withcolor";
        end

        if nPlots == 1
            % legend
            legend(fc + " Hz", 'Orientation', 'vertical', 'FontSize', opts.legendFontSize, Location='eastoutside');
        end
    
        % Image Overlay
        if ~isempty(opts.imgPath)
            axes1 = axes('Parent',fig,'Position',axHor.Position .* [1 1 0.5 0.5] + axHor.Position .* opts.imgPos(1,:));
            axis off
            hold(axes1,'on');
            [img, ~, alphachannel] = imread(opts.imgPath{1});
            alphachannel = alphachannel .* opts.imgAlpha;
            image(img, 'AlphaData', alphachannel);
            axis equal
            box off
        else
            title(opts.subPlotTitles{1})
        end
        plotIdx = plotIdx + 1;
    end

    if opts.plotHorMedFront(2)
        % Median Plane
        subplot(1,nPlots,plotIdx)
        % Polar Plot
        polarplot([elMesh(:,1); pi-flip(elMesh(1:end-1,1))], rmsLevelMedian, 'LineWidth', opts.lineWidth)
        rlim([opts.minPlotLevel(2) 0])
        rticks(opts.minPlotLevel(2):opts.rTickRes(2):0);
        rl = rticklabels;
        rl(:,2) = {'dB'};
        rl = join(rl);
        rticklabels(rl);
        thetalim([0 180])
        axMed = gca;
        axMed.FontSize = opts.axFontSize;
        axMed.Color = [1 1 1 0];
        axMed.ThetaColor = [0 0 0];
        axMed.RColor = [0 0 0];
        axMed.LineStyleOrder = opts.lineStyle;
        if ~isMATLABReleaseOlderThan("R2023a")
            axMed.LineStyleCyclingMethod = "withcolor";
        end

        % Image Overlay
        if ~isempty(opts.imgPath)
            axes2 = axes('Parent',fig,'Position',axMed.Position .* [1 1 0.5 0.5] + axMed.Position .* opts.imgPos(2,:));
            axis off
            hold(axes2,'on');
            [img, ~, alphachannel] = imread(opts.imgPath{2});
            alphachannel = alphachannel .* opts.imgAlpha;
            image(flip(img, 1), 'AlphaData', flip(alphachannel,1));
            axis equal
            box off
        else
            title(opts.subPlotTitles{2})
        end
        plotIdx = plotIdx + 1;
    end

    if opts.plotHorMedFront(3)
        % Frontal Plane
        subplot(1,nPlots,plotIdx)
        % Polar Plot
        polarplot([elMesh(:,1); pi-flip(elMesh(1:end-1,1))], rmsLevelFrontal, 'LineWidth', opts.lineWidth)
        rlim([opts.minPlotLevel(3) 0])
        rticks(opts.minPlotLevel(3):opts.rTickRes(3):0);
        rl = rticklabels;
        rl(:,2) = {'dB'};
        rl = join(rl);
        rticklabels(rl);
        thetalim([0 180])
        axFront = gca;
        axFront.FontSize = opts.axFontSize;
        axFront.Color = [1 1 1 0];
        axFront.ThetaColor = [0 0 0];
        axFront.RColor = [0 0 0];
        axFront.LineStyleOrder = opts.lineStyle;
        if ~isMATLABReleaseOlderThan("R2023a")
            axFront.LineStyleCyclingMethod = "withcolor";
        end

        % legend
        hL = legend(fc + " Hz", 'Orientation', 'horizontal', 'FontSize', opts.legendFontSize);
        hL.Position(1:2) = [axMed.Position(1)*0.95, axMed.Position(2) * 0.8]; 
    
        % Image Overlay
        if ~isempty(opts.imgPath)
            axes3 = axes('Parent',fig,'Position',axFront.Position .* [1 1 0.5 0.5] + axFront.Position .* opts.imgPos(3,:));
            axis off
            hold(axes3,'on');
        
            [img, ~, alphachannel] = imread(opts.imgPath{3});  
            alphachannel = alphachannel .* opts.imgAlpha;
            image(flip(img, 1), 'AlphaData', flip(alphachannel,1));
            axis equal
            box off
        else
            title(opts.subPlotTitles{3})
        end
    end
    
    % Flip order of figure children so that overlay is in the back
    fig.Children = flipud(fig.Children);
end