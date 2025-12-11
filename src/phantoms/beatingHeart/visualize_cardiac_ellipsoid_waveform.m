function visualize_cardiac_ellipsoid_waveform(t_s, a_mm, b_mm, visOpts)
%VISUALIZE_CARDIAC_ELLIPSOID_WAVEFORM  Animate LV ellipsoid semi-axes waveform.
%
%   visualize_cardiac_ellipsoid_waveform(t_s, a_mm, b_mm, visOpts)
%
%   Inputs:
%       t_s     - time [s], 1xN, strictly increasing
%       a_mm    - ellipsoid semi-axis along long-axis [mm], 1xN
%       b_mm    - ellipsoid short-axis [mm], 1xN
%       visOpts - (optional struct)
%                 .frameStep - frame stride for animation (default ~N/500)
%                 .nTheta    - number of angular samples (default 200)
%
%   Visualization:
%       - Coronal ellipse outline (x = b, y = a) animated over time.
%       - Semi-axes a(t) and b(t) vs time with a moving cursor.

    arguments
        t_s    (1,:) double {mustBeReal, mustBeFinite}
        a_mm   (1,:) double {mustBeReal, mustBeFinite, mustBePositive}
        b_mm   (1,:) double {mustBeReal, mustBeFinite, mustBePositive}
        visOpts struct = struct()
    end

    N = numel(t_s);
    if any(diff(t_s) <= 0)
        error('t_s must be strictly increasing.');
    end
    if numel(a_mm) ~= N || numel(b_mm) ~= N
        error('a_mm and b_mm must have the same length as t_s.');
    end

    if ~isfield(visOpts, 'frameStep') || isempty(visOpts.frameStep)
        visOpts.frameStep = max(1, floor(N/500));
    end
    if ~isfield(visOpts, 'nTheta') || isempty(visOpts.nTheta)
        visOpts.nTheta = 200;
    end

    frameStep = visOpts.frameStep;
    nTheta    = visOpts.nTheta;

    a_cm = a_mm / 10;          % convert to cm for plotting
    b_cm = b_mm / 10;

    theta = linspace(0, 2*pi, nTheta);

    maxA = max(a_cm);
    maxB = max(b_cm);
    pad = 1.2;

    figure('Name','Cardiac ellipsoid (simplified)','Color','w');
    tlo = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    % (1) Coronal ellipse outline (x = b, y = a)
    axEll = nexttile(tlo, 1);
    hold(axEll, 'on');
    x0 = b_cm(1) * cos(theta);
    y0 = a_cm(1) * sin(theta);
    hEll = plot(axEll, x0, y0, 'b', 'LineWidth', 2);
    axis(axEll, 'equal');
    xlim(axEll, pad * [-maxB, maxB]);
    ylim(axEll, pad * [-maxA, maxA]);
    xlabel(axEll, 'Left–Right [cm]');
    ylabel(axEll, 'Long-axis [cm]');
    title(axEll, 'Ellipsoid cross-section');
    grid(axEll, 'on');

    % (2) Semi-axes vs time with moving cursor
    axAxes = nexttile(tlo, 2);
    hold(axAxes, 'on');
    plot(axAxes, t_s, a_cm, 'LineWidth', 1.5, 'DisplayName', 'a (long)');
    plot(axAxes, t_s, b_cm, 'LineWidth', 1.5, 'DisplayName', 'b (short)');
    hPointA = plot(axAxes, t_s(1), a_cm(1), 'ok', 'MarkerFaceColor', 'k', ...
        'MarkerSize', 6, 'DisplayName', 'Cursor a');
    hPointB = plot(axAxes, t_s(1), b_cm(1), 'or', 'MarkerFaceColor', 'r', ...
        'MarkerSize', 6, 'DisplayName', 'Cursor b');
    xlabel(axAxes, 'Time [s]');
    ylabel(axAxes, 'Semi-axis [cm]');
    title(axAxes, 'Semi-axes vs time');
    legend(axAxes, 'Location', 'best');
    grid(axAxes, 'on');

    for k = 1:frameStep:numel(t_s)
        ak = a_cm(k);
        bk = b_cm(k);

        xk = bk * cos(theta);
        yk = ak * sin(theta);
        set(hEll, 'XData', xk, 'YData', yk);

        set(hPointA, 'XData', t_s(k), 'YData', ak);
        set(hPointB, 'XData', t_s(k), 'YData', bk);

        drawnow;
    end
end
