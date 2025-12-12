function displayBreathingMotion(t_s, V_L, ellipsoidParams, B_phase, visOpts)
%DISPLAYBREATHINGMOTION  Animate simple lung ellipsoids + circles + V(t).
%
%   displayBreathingMotion(t_s, V_L, ellipsoidParams, B_phase, visOpts)
%
%   Inputs:
%       t_s     - time [s], 1xN
%       V_L     - total lung volume [L], 1xN
%       ellipsoidParams - struct with fields:
%           .axes            - struct('a_mm', R_mm, 'b_mm', R_mm, 'c_mm', H_mm)
%           .leftCenter_mm   - Nx3 centers for left lung
%           .rightCenter_mm  - Nx3 centers for right lung
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
        ellipsoidParams struct
        B_phase (1,:) double {mustBeReal, mustBeFinite}
        visOpts struct = struct()   % <- visOpts is a struct positional arg
    end

    N = numel(t_s);
    if any(diff(t_s) <= 0)
        error('t_s must be strictly increasing.');
    end
    if any([numel(V_L), numel(B_phase)] ~= N)
        error('All input vectors must have the same length as t_s.');
    end

    requiredFields = {'axes', 'leftCenter_mm', 'rightCenter_mm'};
    for idxField = 1:numel(requiredFields)
        if ~isfield(ellipsoidParams, requiredFields{idxField})
            error('visualize_lung_motion:MissingField', ...
                'ellipsoidParams.%s is required.', requiredFields{idxField});
        end
    end

    axesParams = ellipsoidParams.axes;
    for axisField = {'a_mm', 'b_mm', 'c_mm'}
        if ~isfield(axesParams, axisField{1})
            error('visualize_lung_motion:MissingAxisField', ...
                'ellipsoidParams.axes.%s is required.', axisField{1});
        end
    end

    validateattributes(axesParams.a_mm, {'numeric'}, {'real', 'positive'});
    validateattributes(axesParams.b_mm, {'numeric'}, {'real', 'positive'});
    validateattributes(axesParams.c_mm, {'numeric'}, {'real', 'positive'});

    if ~isequal(size(axesParams.a_mm), size(axesParams.b_mm), size(axesParams.c_mm))
        error('visualize_lung_motion:AxisLengthMismatch', ...
            'ellipsoidParams.axes fields must share the same size.');
    end

    R_mm = axesParams.a_mm;
    H_mm = axesParams.c_mm;

    if ~isscalar(R_mm) && numel(R_mm) ~= N
        error('visualize_lung_motion:RadiusLengthMismatch', ...
            'ellipsoidParams.axes values must be scalar or length N.');
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

    leftCenter_mm = ellipsoidParams.leftCenter_mm;
    rightCenter_mm = ellipsoidParams.rightCenter_mm;

    if size(leftCenter_mm, 1) ~= N || size(rightCenter_mm, 1) ~= N
        error('visualize_lung_motion:CenterLengthMismatch', ...
            'Center arrays must have N rows.');
    end

    % Work in cm for display
    R_cm = R_mm / 10;
    H_cm = H_mm / 10;
    leftCenter_cm = leftCenter_mm / 10;
    rightCenter_cm = rightCenter_mm / 10;

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

    xL0 = leftCenter_cm(1,1);
    xR0 = rightCenter_cm(1,1);
    zL0 = leftCenter_cm(1,3);
    zR0 = rightCenter_cm(1,3);

    % Coronal ellipses
    axes(axCor); %#ok<LAXES>
    hold(axCor,'on');
    xL_ell = xL0 + R0 * cos(theta);
    yL_ell = zL0 + (H0/2) * sin(theta);
    hEllL = plot(axCor, xL_ell, yL_ell,'b','LineWidth',2);

    xR_ell = xR0 + R0 * cos(theta);
    yR_ell = zR0 + (H0/2) * sin(theta);
    hEllR = plot(axCor, xR_ell, yR_ell,'r','LineWidth',2);

    axis(axCor,'equal');
    xlabel(axCor,'Left–Right [cm]');
    ylabel(axCor,'Inferior–Superior [cm]');
    title(axCor,'Coronal view (ellipsoids)');

    maxR   = max(R_cm);
    maxH   = max(H_cm);
    maxSep = max(abs([leftCenter_cm(:,1); rightCenter_cm(:,1)])) + maxR;
    maxZ   = max(abs([leftCenter_cm(:,3); rightCenter_cm(:,3)])) + maxH/2;
    xlim(axCor,[-maxSep, maxSep]);
    ylim(axCor,[-maxZ, maxZ]);

    % Axial circles
    axes(axAx); %#ok<LAXES>
    hold(axAx,'on');
    xL_circ = xL0 + R0 * cos(theta);
    yL_circ = leftCenter_cm(1,2) + R0 * sin(theta);
    hCircL = plot(axAx, xL_circ, yL_circ,'b','LineWidth',2);

    xR_circ = xR0 + R0 * cos(theta);
    yR_circ = rightCenter_cm(1,2) + R0 * sin(theta);
    hCircR = plot(axAx, xR_circ, yR_circ,'r','LineWidth',2);

    axis(axAx,'equal');
    xlabel(axAx,'Left–Right [cm]');
    ylabel(axAx,'Anterior–Posterior [cm]');
    title(axAx,'Axial view (circles)');

    maxY = max(abs([leftCenter_cm(:,2); rightCenter_cm(:,2)])) + maxR;
    xlim(axAx,[-maxSep, maxSep]);
    ylim(axAx,[-maxY, maxY]);

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
        xL  = leftCenter_cm(k,1);
        xR  = rightCenter_cm(k,1);
        yL  = leftCenter_cm(k,2);
        yR  = rightCenter_cm(k,2);
        zL  = leftCenter_cm(k,3);
        zR  = rightCenter_cm(k,3);

        % Update coronal ellipses
        xL_ell = xL + Rk * cos(theta);
        yL_ell = zL + (Hk/2) * sin(theta);
        xR_ell = xR + Rk * cos(theta);
        yR_ell = zR + (Hk/2) * sin(theta);
        set(hEllL,'XData',xL_ell,'YData',yL_ell);
        set(hEllR,'XData',xR_ell,'YData',yR_ell);

        % Update axial circles
        xL_circ = xL + Rk * cos(theta);
        yL_circ = yL + Rk * sin(theta);
        xR_circ = xR + Rk * cos(theta);
        yR_circ = yR + Rk * sin(theta);
        set(hCircL,'XData',xL_circ,'YData',yL_circ);
        set(hCircR,'XData',xR_circ,'YData',yR_circ);

        % Update volume cursor
        set(tLineV,'XData',[tk tk],'YData',ylV);

        drawnow;
    end
end
