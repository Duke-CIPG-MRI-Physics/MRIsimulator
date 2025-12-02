function visualize_cardiac_ellipsoid(t_s, V_ml, a_m, b_m, c_m, eps_L, eps_C, visOpts)
%VISUALIZE_CARDIAC_ELLIPSOID  Animate LV ellipsoid + volume + strain.
%
%   visualize_cardiac_ellipsoid(t_s, V_ml, a_m, b_m, c_m, eps_L, eps_C, visOpts)
%
%   Inputs:
%       t_s     - time [s], 1xN, strictly increasing
%       V_ml    - LV volume [mL] vs time, 1xN
%       a_m     - ellipsoid semi-axis along long-axis [m], 1xN
%       b_m     - ellipsoid short-axis [m], 1xN
%       c_m     - ellipsoid short-axis [m], 1xN (not used separately here)
%       eps_L   - longitudinal strain (GLS, relative to ED), 1xN
%       eps_C   - circumferential strain (GCS, relative to ED), 1xN
%
%       visOpts (struct, optional):
%           .frameStep - subsampling factor for animation (default ~300 frames)
%           .nTheta    - number of angular samples for ellipse/circle (default 200)
%
%   Layout:
%       Subplot (2,2,1): Coronal view   (ellipse: x=b, y=a)
%       Subplot (2,2,2): Axial view     (circle:  radius=b)
%       Subplot (2,2,3): LV volume vs time   + moving point cursor
%       Subplot (2,2,4): Global strains vs time (GLS, GCS) + moving point cursors

    arguments
        t_s    (1,:) double {mustBeReal, mustBeFinite}
        V_ml   (1,:) double {mustBeReal, mustBeFinite}
        a_m    (1,:) double {mustBeReal, mustBeFinite, mustBePositive}
        b_m    (1,:) double {mustBeReal, mustBeFinite, mustBePositive}
        c_m    (1,:) double {mustBeReal, mustBeFinite, mustBePositive} %#ok<INUSD>
        eps_L  (1,:) double {mustBeReal, mustBeFinite}
        eps_C  (1,:) double {mustBeReal, mustBeFinite}
        visOpts struct = struct()
    end

    % ---------------------------------------------------------------------
    % Defaults for visOpts
    % ---------------------------------------------------------------------
    N = numel(t_s);
    if ~isfield(visOpts,'frameStep') || isempty(visOpts.frameStep)
        visOpts.frameStep = max(1, floor(N/300));   % ~300 frames max
    end
    if ~isfield(visOpts,'nTheta') || isempty(visOpts.nTheta)
        visOpts.nTheta = 200;
    end
    frameStep = visOpts.frameStep;
    nTheta    = visOpts.nTheta;

    % Basic checks
    if any(diff(t_s) <= 0)
        error('t_s must be strictly increasing.');
    end
    if any([numel(V_ml), numel(a_m), numel(b_m), numel(eps_L), numel(eps_C)] ~= N)
        error('All input vectors must have the same length as t_s.');
    end

    % ---------------------------------------------------------------------
    % Pre-compute geometry in cm for plotting
    % ---------------------------------------------------------------------
    a_cm = 100 * a_m;
    b_cm = 100 * b_m;

    theta = linspace(0, 2*pi, nTheta);

    % Axis limits
    maxA = max(a_cm);
    maxB = max(b_cm);
    Rmax = maxB;    % axial radius
    padA = 1.2;
    padB = 1.2;

    % ---------------------------------------------------------------------
    % Figure + tiled layout 2x2
    % ---------------------------------------------------------------------
    figure('Name','Cardiac Ellipsoid Visualization','Color','w');
    tlo = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    axCor   = nexttile(tlo,1);   % (2,2,1) coronal
    axAxial = nexttile(tlo,2);   % (2,2,2) axial
    axV     = nexttile(tlo,3);   % (2,2,3) volume
    axStr   = nexttile(tlo,4);   % (2,2,4) strains

    % ---------------------------------------------------------------------
    % (1) Coronal view: ellipse (x=b, y=a)
    % ---------------------------------------------------------------------
    axes(axCor); %#ok<LAXES>
    hold(axCor,'on');
    xCor0 = b_cm(1) * cos(theta);
    yCor0 = a_cm(1) * sin(theta);
    hEllCor = plot(axCor, xCor0, yCor0, 'b', 'LineWidth', 2);
    axis(axCor,'equal');
    xlabel(axCor,'Left–Right [cm]');
    ylabel(axCor,'Long-axis [cm]');
    title(axCor,'Coronal view (ellipsoid)');
    xlim(axCor, padB * [-maxB, maxB]);
    ylim(axCor, padA * [-maxA, maxA]);
    grid(axCor,'on');

    % ---------------------------------------------------------------------
    % (2) Axial view: circle (short-axis vs short-axis, radius=b)
    % ---------------------------------------------------------------------
    axes(axAxial); %#ok<LAXES>
    hold(axAxial,'on');
    xAx0 = b_cm(1) * cos(theta);
    yAx0 = b_cm(1) * sin(theta);
    hCircAx = plot(axAxial, xAx0, yAx0, 'r', 'LineWidth', 2);
    axis(axAxial,'equal');
    xlabel(axAxial,'Short-axis X [cm]');
    ylabel(axAxial,'Short-axis Y [cm]');
    title(axAxial,'Axial view (short-axis circle)');
    xlim(axAxial, padB * [-Rmax, Rmax]);
    ylim(axAxial, padB * [-Rmax, Rmax]);
    grid(axAxial,'on');

    % ---------------------------------------------------------------------
    % (3) LV volume vs time with moving point cursor
    % ---------------------------------------------------------------------
    axes(axV);
    hold(axV,'on');
    plot(axV, t_s, V_ml, 'LineWidth', 1.5);
    xlabel(axV,'Time [s]');
    ylabel(axV,'V_{LV} [mL]');
    title(axV,'LV volume vs time');
    grid(axV,'on');
    ylV = ylim(axV); %#ok<NASGU>  % keep for consistency if you tweak later

    % Point cursor for volume
    hPointV = plot(axV, t_s(1), V_ml(1), 'ok', ...
                   'MarkerFaceColor','k','MarkerSize',6);

    % ---------------------------------------------------------------------
    % (4) Global strains vs time with moving point cursors
    % ---------------------------------------------------------------------
    axes(axStr);
    hold(axStr,'on');
    hStrL = plot(axStr, t_s, eps_L*100, 'LineWidth',1.5); hold(axStr,'on');
    hStrC = plot(axStr, t_s, eps_C*100, 'LineWidth',1.5);
    xlabel(axStr,'Time [s]');
    ylabel(axStr,'Strain [%]');
    legend(axStr, {'GLS','GCS'}, 'Location','best');
    title(axStr,'Global strains vs time');
    grid(axStr,'on');
    ylStr = ylim(axStr); %#ok<NASGU>

    % Point cursors for GLS and GCS
    hPointStrL = plot(axStr, t_s(1), eps_L(1)*100, 'ok', ...
                      'MarkerFaceColor','k','MarkerSize',6);
    hPointStrC = plot(axStr, t_s(1), eps_C(1)*100, 'or', ...
                      'MarkerFaceColor','r','MarkerSize',6);

    % Make sure lines are behind the point markers
    uistack(hStrL,'bottom');
    uistack(hStrC,'bottom');

    % ---------------------------------------------------------------------
    % Animation loop
    % ---------------------------------------------------------------------
    for k = 1:frameStep:N
        tk  = t_s(k);
        ak  = a_cm(k);
        bk  = b_cm(k);

        % Update coronal ellipse
        xCor = bk * cos(theta);
        yCor = ak * sin(theta);
        set(hEllCor, 'XData', xCor, 'YData', yCor);

        % Update axial circle
        xAx = bk * cos(theta);
        yAx = bk * sin(theta);
        set(hCircAx, 'XData', xAx, 'YData', yAx);

        % Update volume point cursor
        set(hPointV, 'XData', tk, 'YData', V_ml(k));

        % Update strain point cursors
        set(hPointStrL, 'XData', tk, 'YData', eps_L(k)*100);
        set(hPointStrC, 'XData', tk, 'YData', eps_C(k)*100);

        drawnow;
    end
end
