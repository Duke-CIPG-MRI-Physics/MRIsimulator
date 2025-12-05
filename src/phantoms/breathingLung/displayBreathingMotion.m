function displayBreathingMotion(t_s, V_L, R_mm, H_mm, B_phase, visOpts)
%DISPLAYBREATHINGMOTION  Animate simple lung ellipsoids + circles + V(t).
%
%   displayBreathingMotion(t_s, V_L, R_mm, H_mm, B_phase, visOpts)
%
%   Inputs:
%       t_s     - time [s], 1xN
%       V_L     - total lung volume [L], 1xN
%       R_mm    - lung radius [mm], 1xN
%       H_mm    - lung height [mm], 1xN
%       B_phase - breathing phase [rad], 1xN (for consistency; not used)
%       visOpts - (optional struct)
%                 .frameStep - frame stride for animation
%                 .nTheta    - number of points for ellipse/circle
%
%   Visualization:
%       - Coronal: two ellipses per frame (width ~ 2*R, height = H)
%       - Axial:   two circles per frame (radius = R)
%       - V(t) panel with a moving cursor.

    arguments
        t_s     (1,:) double {mustBeReal, mustBeFinite}
        V_L     (1,:) double {mustBeReal, mustBeFinite}
        R_mm    (1,:) double {mustBeReal, mustBeFinite, mustBePositive}
        H_mm    (1,:) double {mustBeReal, mustBeFinite, mustBePositive}
        B_phase (1,:) double {mustBeReal, mustBeFinite}
        visOpts struct = struct()   % <- visOpts is a struct positional arg
    end

    N = numel(t_s);
    if any(diff(t_s) <= 0)
        error('t_s must be strictly increasing.');
    end
    if any([numel(V_L), numel(R_mm), numel(H_mm), numel(B_phase)] ~= N)
        error('All input vectors must have the same length as t_s.');
    end

    % Defaults for visualization options
    if ~isfield(visOpts,'frameStep') || isempty(visOpts.frameStep)
        visOpts.frameStep = max(1, floor(N/300));   % ~300 frames max
    end
    if ~isfield(visOpts,'nTheta') || isempty(visOpts.nTheta)
        visOpts.nTheta = 200;
    end

    frameStep = visOpts.frameStep;
    nTheta    = visOpts.nTheta;

    % Work in cm for display
    R_cm = R_mm / 10;
    H_cm = H_mm / 10;

    % Angular sampling for circles/ellipses
    theta = linspace(0, 2*pi, nTheta);

    % ---------------------------------------------------------------------
    % Figure + tiled layout
    % ---------------------------------------------------------------------
    figure('Name','Breathing Motion Visualization','Color','w');
    % 2x2 layout: [Coronal, Axial; V(t), (unused)]
    tlo = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    axCor = nexttile(tlo,1);
    axAx  = nexttile(tlo,2);
    axV   = nexttile(tlo,3);
    % tile 4 unused for now

    % ---------------------------------------------------------------------
    % Initial geometry
    % ---------------------------------------------------------------------
    R0   = R_cm(1);
    H0   = H_cm(1);
    sep0 = 3 * R0;     % center-to-center separation in cm

    xL0 = -sep0/2;
    xR0 =  sep0/2;

    % Coronal ellipses
    axes(axCor); %#ok<LAXES>
    hold(axCor,'on');
    xL_ell = xL0 + R0 * cos(theta);
    yL_ell =       (H0/2) * sin(theta);
    hEllL = plot(axCor, xL_ell, yL_ell,'b','LineWidth',2);

    xR_ell = xR0 + R0 * cos(theta);
    yR_ell =       (H0/2) * sin(theta);
    hEllR = plot(axCor, xR_ell, yR_ell,'r','LineWidth',2);

    axis(axCor,'equal');
    xlabel(axCor,'Left–Right [cm]');
    ylabel(axCor,'Inferior–Superior [cm]');
    title(axCor,'Coronal view (ellipsoids)');

    maxR   = max(R_cm);
    maxH   = max(H_cm);
    maxSep = 3 * maxR;
    xlim(axCor,[-maxSep, maxSep]);
    ylim(axCor,[-maxSep, maxSep]);

    % Axial circles
    axes(axAx); %#ok<LAXES>
    hold(axAx,'on');
    xL_circ = xL0 + R0 * cos(theta);
    yL_circ =       R0 * sin(theta);
    hCircL = plot(axAx, xL_circ, yL_circ,'b','LineWidth',2);

    xR_circ = xR0 + R0 * cos(theta);
    yR_circ =       R0 * sin(theta);
    hCircR = plot(axAx, xR_circ, yR_circ,'r','LineWidth',2);

    axis(axAx,'equal');
    xlabel(axAx,'Left–Right [cm]');
    ylabel(axAx,'Anterior–Posterior [cm]');
    title(axAx,'Axial view (circles)');

    xlim(axAx,[-maxSep, maxSep]);
    ylim(axAx,[-maxSep, maxSep]);

    % Volume vs time
    axes(axV);
    hold(axV,'on');
    plot(axV, t_s, V_L,'LineWidth',1.5);
    ylabel(axV,'V_{lung} [L]');
    xlabel(axV,'Time [s]');
    title(axV,'Lung volume vs time');
    grid(axV,'on');
    ylV = ylim(axV);
    tLineV = line(axV,[t_s(1) t_s(1)], ylV, ...
        'LineStyle','--','Color',[0 0 0],'LineWidth',1);

    % ---------------------------------------------------------------------
    % Animation loop
    % ---------------------------------------------------------------------
    for k = 1:frameStep:N
        tk  = t_s(k);
        Rk  = R_cm(k);
        Hk  = H_cm(k);
        sep = 3 * Rk;
        xL  = -sep/2;
        xR  =  sep/2;

        % Update coronal ellipses
        xL_ell = xL + Rk * cos(theta);
        yL_ell =      (Hk/2) * sin(theta);
        xR_ell = xR + Rk * cos(theta);
        yR_ell =      (Hk/2) * sin(theta);
        set(hEllL,'XData',xL_ell,'YData',yL_ell);
        set(hEllR,'XData',xR_ell,'YData',yR_ell);

        % Update axial circles
        xL_circ = xL + Rk * cos(theta);
        yL_circ =      Rk * sin(theta);
        xR_circ = xR + Rk * cos(theta);
        yR_circ =      Rk * sin(theta);
        set(hCircL,'XData',xL_circ,'YData',yL_circ);
        set(hCircR,'XData',xR_circ,'YData',yR_circ);

        % Update volume cursor
        set(tLineV,'XData',[tk tk],'YData',ylV);

        drawnow;
    end
end
