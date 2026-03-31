function plotTWISTFrameOrdering( ...
    Sampling_Table, Matrix_Size_Acquired, regionA, plotTitle, maxFramesToShow)
% plotTWISTFrameOrdering  Visualize TWIST B-frame ordering on the PE plane.
%   plotTWISTFrameOrdering(Sampling_Table, Matrix_Size_Acquired, regionA,
%       plotTitle, maxFramesToShow) creates a tiled figure showing the
%   acquisition trajectory for the partial B frames. The trajectory is drawn
%   on the phase/slice plane, with region A shown in gray and B samples
%   colored by in-frame acquisition order.

    arguments
        Sampling_Table table
        Matrix_Size_Acquired (1,3) {mustBeNumeric, mustBePositive, mustBeInteger}
        regionA = []
        plotTitle (1,1) string = "TWIST Frame Ordering"
        maxFramesToShow (1,1) {mustBeNumeric, mustBePositive, mustBeInteger} = 9
    end

    if ismember("Frequency", Sampling_Table.Properties.VariableNames)
        planeTable = Sampling_Table(Sampling_Table.Frequency == 1, :);
    else
        planeTable = Sampling_Table;
    end

    if isempty(planeTable)
        warning("TWIST:plotTWISTFrameOrdering:EmptyTable", ...
            "Sampling table is empty. Nothing to visualize.");
        return
    end

    frames = unique(planeTable.Frame, "stable")';
    if ismember("Bj", planeTable.Properties.VariableNames)
        nBSubsets = max(planeTable.Bj);
    else
        nBSubsets = 0;
    end

    partialFrames = selectPartialFrames(planeTable, frames, nBSubsets);
    framesToShow = partialFrames(1:min(maxFramesToShow, numel(partialFrames)));

    nTiles = ceil(sqrt(numel(framesToShow)));
    figure('Name', char(plotTitle), 'Color', 'w');
    tiledlayout(nTiles, nTiles, 'TileSpacing', 'compact', 'Padding', 'compact');

    [aRows, aCols] = findRegionAMarkers(regionA, planeTable);
    for iFrame = 1:numel(framesToShow)
        currentFrame = framesToShow(iFrame);
        frameRows = planeTable(planeTable.Frame == currentFrame, :);
        bRows = frameRows(frameRows.Bj > 0, :);

        nexttile;
        hold on;
        scatter(aRows, aCols, 10, repmat(0.85, numel(aRows), 3), 'filled');

        if ~isempty(bRows)
            xData = bRows.("Row (phase)");
            yData = bRows.("Column (slice)");
            plot(xData, yData, '-', 'Color', [0.2, 0.2, 0.2], 'LineWidth', 0.5);
            scatter(xData, yData, 22, 1:height(bRows), 'filled');
            plot(xData(1), yData(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 6);
            plot(xData(end), yData(end), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
        end

        axis image;
        xlim([1, Matrix_Size_Acquired(2)]);
        ylim([1, Matrix_Size_Acquired(3)]);
        set(gca, 'YDir', 'normal');
        title(buildFrameTitle(currentFrame, bRows), 'Interpreter', 'none');
        xlabel('Phase Encode');
        ylabel('Slice Encode');
        grid on;
        box on;
    end

    colormap(turbo);
    sgtitle(plotTitle);
end

function partialFrames = selectPartialFrames(planeTable, frames, nBSubsets)
% selectPartialFrames  Prefer partial B frames over full anchor frames.

    partialFrames = [];
    for iFrame = 1:numel(frames)
        frameRows = planeTable(planeTable.Frame == frames(iFrame), :);
        if nBSubsets == 0
            partialFrames(end + 1) = frames(iFrame); %#ok<AGROW>
            continue
        end

        presentSubsets = unique(frameRows.Bj(frameRows.Bj > 0))';
        if numel(presentSubsets) < nBSubsets
            partialFrames(end + 1) = frames(iFrame); %#ok<AGROW>
        end
    end

    if isempty(partialFrames)
        partialFrames = frames;
    end
end

function [aRows, aCols] = findRegionAMarkers(regionA, planeTable)
% findRegionAMarkers  Determine region-A markers for the plot background.

    if ~isempty(regionA)
        [aRows, aCols] = find(regionA);
    else
        aRows = planeTable.("Row (phase)")(planeTable.Bj == 0);
        aCols = planeTable.("Column (slice)")(planeTable.Bj == 0);
    end
end

function titleText = buildFrameTitle(frameNumber, bRows)
% buildFrameTitle  Describe one frame's active Bj subset(s).

    if isempty(bRows)
        titleText = sprintf('Frame %d | A only', frameNumber);
        return
    end

    bjLabels = unique(bRows.Bj)';
    titleText = sprintf('Frame %d | B%s', frameNumber, strjoin(string(bjLabels), ','));
end
