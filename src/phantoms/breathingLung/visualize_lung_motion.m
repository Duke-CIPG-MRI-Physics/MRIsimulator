function displayBreathingMotion(t_s, leftEllipsoidParams, rightEllipsoidParams, visOpts)
%DISPLAYBREATHINGMOTION  Animate simple lung ellipsoids + circles + V(t).
%
%   displayBreathingMotion(t_s, leftEllipsoidParams, rightEllipsoidParams, visOpts)
%
%   Inputs:
%       t_s     - time [s], 1xN
%       leftEllipsoidParams, rightEllipsoidParams - structs with fields:
%           .a_mm, .b_mm, .c_mm - semi-axis lengths [mm]
%           .pose              - pose struct including center.x_mm, etc.
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
        leftEllipsoidParams struct
        rightEllipsoidParams struct
        visOpts struct = struct()   % <- visOpts is a struct positional arg
    end

    N = numel(t_s);
    if any(diff(t_s) <= 0)
        error('t_s must be strictly increasing.');
    end

    leftRequired = {'a_mm', 'b_mm', 'c_mm'};
    rightRequired = {'a_mm', 'b_mm', 'c_mm'};

    for idxField = 1:numel(leftRequired)
        if ~isfield(leftEllipsoidParams, leftRequired{idxField})
            error('visualize_lung_motion:MissingField', ...
                'leftEllipsoidParams.%s is required.', leftRequired{idxField});
        end
    end

    for idxField = 1:numel(rightRequired)
        if ~isfield(rightEllipsoidParams, rightRequired{idxField})
            error('visualize_lung_motion:MissingField', ...
                'rightEllipsoidParams.%s is required.', rightRequired{idxField});
        end
    end

    requiredPoseFields = {'pose'};
    for idxField = 1:numel(requiredPoseFields)
        if ~isfield(leftEllipsoidParams, requiredPoseFields{idxField}) || ...
                ~isstruct(leftEllipsoidParams.pose)
            error('visualize_lung_motion:MissingPose', ...
                'leftEllipsoidParams.pose is required.');
        end
        if ~isfield(rightEllipsoidParams, requiredPoseFields{idxField}) || ...
                ~isstruct(rightEllipsoidParams.pose)
            error('visualize_lung_motion:MissingPose', ...
                'rightEllipsoidParams.pose is required.');
        end
    end

    centerFields = {'x_mm', 'y_mm', 'z_mm'};
    for idxField = 1:numel(centerFields)
        if ~isfield(leftEllipsoidParams.pose, 'center') || ...
                ~isstruct(leftEllipsoidParams.pose.center) || ...
                ~isfield(leftEllipsoidParams.pose.center, centerFields{idxField})
            error('visualize_lung_motion:MissingCenter', ...
                'leftEllipsoidParams.pose.center.%s is required.', centerFields{idxField});
        end

        if ~isfield(rightEllipsoidParams.pose, 'center') || ...
                ~isstruct(rightEllipsoidParams.pose.center) || ...
                ~isfield(rightEllipsoidParams.pose.center, centerFields{idxField})
            error('visualize_lung_motion:MissingCenter', ...
                'rightEllipsoidParams.pose.center.%s is required.', centerFields{idxField});
        end
    end

    R_mm = leftEllipsoidParams.a_mm;
    H_mm = leftEllipsoidParams.c_mm;

    validateattributes(R_mm, {'numeric'}, {'real', 'positive'});
    validateattributes(H_mm, {'numeric'}, {'real', 'positive'});

    if ~isequal(size(R_mm), size(H_mm))
        error('visualize_lung_motion:AxisLengthMismatch', ...
            'ellipsoidParams fields must share the same size.');
    end

    if ~isscalar(R_mm) && numel(R_mm) ~= N
        error('visualize_lung_motion:RadiusLengthMismatch', ...
            'ellipsoidParams values must be scalar or length N.');
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

    leftCenter_mm = [leftEllipsoidParams.pose.center.x_mm(:), ...
        leftEllipsoidParams.pose.center.y_mm(:), ...
        leftEllipsoidParams.pose.center.z_mm(:)];
    rightCenter_mm = [rightEllipsoidParams.pose.center.x_mm(:), ...
        rightEllipsoidParams.pose.center.y_mm(:), ...
        rightEllipsoidParams.pose.center.z_mm(:)];

    % Work in cm for display
    R_cm = R_mm / 10;
    H_cm = H_mm / 10;
    leftCenter_cm = leftCenter_mm / 10;
    rightCenter_cm = rightCenter_mm / 10;

    % Compute volume from geometry for display (both lungs combined)
    V_L = (4/3) * pi .* R_cm.^2 .* H_cm / 1000;

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
