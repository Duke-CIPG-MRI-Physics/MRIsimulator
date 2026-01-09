clc; clear; close all;

%% --- Parameters and Initialization ---
kspaceSize = [256 256 256];  % [frequency, phase (rows), slice (columns)]
pA = 0.016;                  % Proportion of A region (central k-space)
N = 10;                     % Number of time frames (B interleaves)
TR = 6e-3; %not critical, just gives us temporal resolution
%% --- Coordinate Grid Setup (Phase/Slice) ---
kyi = 1:kspaceSize(3);  % slice (columns)
kzi = 1:kspaceSize(2);  % phase (rows)

centerPixel = ceil(0.5 * (kspaceSize + 1));  % Use full kspaceSize as per your request
ky = kyi - centerPixel(3);  % col offset (slice)
kz = kzi - centerPixel(2);  % row offset (phase)
[kyM, kzM] = meshgrid(ky, kz);  % rows, cols

[theta, kr] = cart2pol(kyM, kzM);
%Fix theta range, direction, and 0
theta(theta>0) = 2*pi-(theta(theta>0)); 
theta = abs(theta);

%kr = round(kr/1)*1;  % Integer bins for radius

%% --- Define A and B Sampling Regions ---
max_kz_ky = max(kr(centerPixel(2), end), kr(end, centerPixel(3)));
kc = pA * max_kz_ky;

regionA = (kr < kc);
regionB = ~regionA;

figure; imshow(regionA, []);
title('Regions "A" (central) and "B" (peripheral)');

%% --- Radial and Angular Sorting ---
columnNames = {'Index','Kr','Theta','Region A?'};
unsortedData = table((1:numel(kr))', kr(:), theta(:), regionA(:),'VariableNames',columnNames);
sortedData = sortrows(unsortedData, [2 3], 'ascend');


%% --- Add Bj Column to Sorted Data ---
% Initialize the 'Bj' column with zeros for all rows.
sortedData.Bj = zeros(height(sortedData), 1);

% Identify the rows that belong to Region B (where 'Region A?' is false).
% These are the rows where 'Bj' needs to count up sequentially.
regionB_rows = find(~sortedData.('Region A?'));

% Calculate the total number of points in Region B.
numRegionB_points = length(regionB_rows);

% Create a repeating sequence from 1 to N for Region B points.
% The 'mod' function, combined with adding 1, creates a sequence that
% cycles from 1 to N.

bj_sequence = mod(0:(numRegionB_points - 1), N) + 1;

% Assign the generated 'bj_sequence' to the 'Bj' column for the
% identified Region B rows.
sortedData.Bj(regionB_rows) = bj_sequence';

sortedData_regionA = sortedData(sortedData.("Region A?"),:);
sortedData_regionB = sortedData(~sortedData.("Region A?"),:);


%% --- Displaying the Frames

% Display full acquistion
filled_kSpace = zeros(kspaceSize(2:3));
filled_kSpace(sortedData_regionB.Index) = sortedData_regionB.Bj;

figure
imagesc(filled_kSpace)
colormap('copper')
title('TWIST Acquisition Frame Assignment');

%Display each frame as layers
kspace_full_sequence = zeros([kspaceSize(2:3),N]);
for ii = 1:N
    kspace_full_sequence(:,:,ii) = (filled_kSpace == ii);
end

figure
imshow3D(kspace_full_sequence)
title('Real Frames')
%% --- Creating actual sampling order

%Sampling of k-space starts at the outer edge of A (kr≈Kc) and proceeds 
% toward the origin...
sortedData_regionA_descend = flipud(sortedData_regionA);

%...via all the odd (kr,Θ) points from the sorted list.
A_to_origin = sortedData_regionA_descend(1:2:end,:);

%Upon reaching the minimum kr, the sampling direction is reversed 
% and every even point is acquired until the edge of A is reached.
A_outward = flipud(sortedData_regionA_descend(2:2:end,:));

%Start building a new table which is sorted by sampling order
kspaceSamplingOrder_A = [A_to_origin;A_outward];
kspaceSamplingOrder_frames = [];
kspaceSamplingOrder_B = [];

for ii = 1:N
    B_current_Bj = sortedData_regionB(sortedData_regionB.Bj == ii,:);

    %Then in region B the trajectory Bj is sampled by first acquiring 
% every odd point on the way toward larger kr...
    B_outward_current_Bj = B_current_Bj(1:2:end,:);
   
    %and then every even point on the way back toward Kc 
    B_inward_current_Bj = flipud(B_current_Bj(2:2:end,:));

    kspaceSamplingOrder_B_current_Bj = [B_outward_current_Bj ; B_inward_current_Bj];

    kspaceSamplingOrder_A.Bj = kspaceSamplingOrder_A.Bj + 1;

    kspaceSamplingOrder_B = [kspaceSamplingOrder_B;kspaceSamplingOrder_B_current_Bj];

    kspaceSamplingOrder_frames = [kspaceSamplingOrder_frames;kspaceSamplingOrder_A;kspaceSamplingOrder_B_current_Bj];
    
end

filled_kSpace = zeros(kspaceSize(2:3));
filled_kSpace(kspaceSamplingOrder_A.Index) = 1:length(kspaceSamplingOrder_A.Index);

figure
imshow(filled_kSpace,[])
title('TWIST Region A Sampling Order');

%% -- Appending the initial full k-space acquisition

%assuming the overall trajectory is the same, but we only collect A at the
%very end

kspaceSamplingOrder_initial = [kspaceSamplingOrder_B;kspaceSamplingOrder_A];
kspaceSamplingOrder_initial.Bj = zeros(height(kspaceSamplingOrder_initial),1);

kspaceSamplingOrder_full = [kspaceSamplingOrder_initial;kspaceSamplingOrder_frames];

kspaceSamplingOrder_full.Order = (1:height(kspaceSamplingOrder_full))';

%% -- Adding columns for row,column indices

[row,col] = ind2sub(kspaceSize(2:3),kspaceSamplingOrder_full.("Index"));
kspaceSamplingOrder_full.("Row (phase)") = row;
kspaceSamplingOrder_full.("Column (slice)") = col;


%% --- TWIST Time
%Just following the paper here
TA_full = kspaceSize(2)*kspaceSize(3)*TR;
S_a = sum(regionA(:));
S_pe = numel(regionA(:));
pB = 1/N;

TA_twist =  TA_full*((S_a/S_pe)+((1-(S_a/S_pe))*pB))

Wait_Time = kspaceSize(2)*kspaceSize(3)*TR;
fprintf('The initial acquisition (wait time) will take %g seconds before TWIST sampling begins\n',Wait_Time)

Temporal_Resolution = (height(kspaceSamplingOrder_A)+height(kspaceSamplingOrder_B_current_Bj))*TR;
fprintf('The temporal resolution of this sequence is %g seconds, not including other techniques\n',Temporal_Resolution)


%% --- Animate the Final Sampling Trajectory (Adjusted Colormap) ---

%Pick which one to animate
animation_table = kspaceSamplingOrder_full;


% This animation visualizes the complete, ordered acquisition sequence.
% It uses a custom colormap where Region A (negative Bj values) is mapped
% to 'lines' and Region B (positive Bj values) is mapped to 'lines'.
% Unfilled pixels are displayed as black.

% --- 1. Animation Setup ---
% Create a new figure for the animation
figure;
set(gcf, 'Position', [150, 150, 700, 550]); % [left, bottom, width, height]

% Determine the full range of Bj values for the colormap
min_Bj = min(animation_table.Bj); % Will be negative
max_Bj = max(animation_table.Bj); % Will be N = 10

% --- 2. Create a Custom Colormap (User Specified) ---
% We need a colormap that maps negative values (Region A) to lines
% and positive values (Region B) to lines.
% Create 'lines' colors for the Region A passes (negative Bj values)
num_A_passes = abs(min_Bj);
cmap_A = lines(num_A_passes);

% Create 'lines' colors for the N Region B frames (positive Bj values)
% Note: The original code used 'lines' for cmap_B, it has been changed to 'lines'
% as per the comment description.
cmap_B = flipud(lines(max_Bj));

% Combine the colormaps. The order is important. cmap_A corresponds to the
% lowest (negative) values, and cmap_B corresponds to the highest (positive) values.
custom_colormap = [cmap_A;0,0,0; cmap_B];

% Apply the colormap to the figure
colormap(custom_colormap);
colormap("lines")

% --- 3. Prepare the Plot Area ---
% Initialize an empty k-space grid with NaN (Not-a-Number).
kspace_animation = nan(kspaceSize(2:3));

% Display the initial empty grid and get a handle to the image object.
h_img = imagesc(kspace_animation);

% *** SET BLACK BACKGROUND AND TRANSPARENCY FOR UNFILLED PIXELS ***
% Get the handle to the current axes and set its color to black ('k').
% This will be the color of the unfilled pixels.
ax = gca;
ax.Color = 'k';

% Set the AlphaData to make NaN values transparent.
% 'AlphaData' is a matrix of the same size as 'CData' (our kspace_animation).
% A value of 1 is opaque, and 0 is transparent.
% By setting AlphaData to be the logical negation of isnan(kspace_animation),
% all NaN values will have an alpha of 0 (transparent) and all other values
% will have an alpha of 1 (opaque).
set(h_img, 'AlphaData', ~isnan(kspace_animation));

% Set plot properties
axis equal tight;
title('k-space Sampling Animation (Acquisition Point: 0)');
xlabel('Slice Encoding (ky)');
ylabel('Phase Encoding (kz)');

% Setup the colorbar to correctly display the frame associations
h_bar = colorbar;
ylabel(h_bar, 'Time Frame (Region B) / Pass (Region A)');
clim([min_Bj, max_Bj]); % Set color limits to the full data range

% --- 4. Animation Loop ---
% Loop through every single point in the final, ordered sampling table.
num_points_total = height(animation_table);

for t = 1:num_points_total

% Get the linear index of the current point to sample
current_index = animation_table.Index(t);

% Get the Bj value, which determines the point's color
current_bj = animation_table.Bj(t);

% Place the point in our animation grid.
kspace_animation(current_index) = current_bj;

% Update the 'CData' of the existing image plot for efficiency.
set(h_img, 'CData', kspace_animation);

% *** UPDATE ALPHADATA IN THE LOOP ***
% As we fill in the kspace_animation matrix, we also need to update
% the AlphaData to make the newly filled pixels opaque.
set(h_img, 'AlphaData', ~isnan(kspace_animation));

% Update the title to show the progress of the acquisition
title('k-space Sampling Animation');

% Force MATLAB to draw the update.
drawnow limitrate;

pause(.0001)

    if t>1 && t<num_points_total
       if kspaceSamplingOrder_full.Bj(t) == (kspaceSamplingOrder_full.Bj(t-1)+1)
            pause(1)
       end
    end
end

title('k-space Sampling Animation - Finished');

% %% --- Animate the Final Sampling Trajectory (and Record Video) ---
% % This animation visualizes the complete, ordered acquisition sequence.
% % It uses a custom colormap where Region A (negative Bj values) is mapped
% % to 'lines' and Region B (positive Bj values) is mapped to 'lines'.
% % Unfilled pixels are displayed as black.
% % The output is recorded to a video file named 'kspace_animation.mp4'.
% 
% % --- 1. Animation and Video Setup ---
% % Create a new figure for the animation
% figure;
% set(gcf, 'Position', [150, 150, 700, 550]); % [left, bottom, width, height]
% 
% % *** VIDEO RECORDING SETUP ***
% % Create and configure the VideoWriter object
% videoFilename = 'kspace_animation.mp4';
% v = VideoWriter(videoFilename, 'MPEG-4');
% v.FrameRate = 60;  % Adjust frame rate for desired speed (e.g., 30, 60)
% open(v); % Open the file for writing
% 
% % Determine the full range of Bj values for the colormap
% min_Bj = min(animation_table.Bj); % Will be negative
% max_Bj = max(animation_table.Bj); % Will be N = 10
% 
% % --- 2. Create a Custom Colormap (User Specified) ---
% % We need a colormap that maps negative values (Region A) to lines
% % and positive values (Region B) to lines.
% % Create 'lines' colors for the Region A passes (negative Bj values)
% num_A_passes = abs(min_Bj);
% cmap_A = lines(num_A_passes);
% 
% % Create 'lines' colors for the N Region B frames (positive Bj values)
% % Note: The original code used 'lines' for cmap_B, it has been changed to 'lines'
% % as per the comment description.
% cmap_B = flipud(lines(max_Bj));
% 
% % Combine the colormaps. The order is important. cmap_A corresponds to the
% % lowest (negative) values, and cmap_B corresponds to the highest (positive) values.
% custom_colormap = [cmap_A;0,0,0; cmap_B];
% 
% % Apply the colormap to the figure
% colormap(custom_colormap);
% colormap("lines")
% 
% % --- 3. Prepare the Plot Area ---
% % Initialize an empty k-space grid with NaN (Not-a-Number).
% kspace_animation = nan(kspaceSize(2:3));
% 
% % Display the initial empty grid and get a handle to the image object.
% h_img = imagesc(kspace_animation);
% 
% % *** SET BLACK BACKGROUND AND TRANSPARENCY FOR UNFILLED PIXELS ***
% % Get the handle to the current axes and set its color to black ('k').
% % This will be the color of the unfilled pixels.
% ax = gca;
% ax.Color = 'k';
% 
% % Set the AlphaData to make NaN values transparent.
% set(h_img, 'AlphaData', ~isnan(kspace_animation));
% 
% % Set plot properties
% axis equal tight;
% title('k-space Sampling Animation');
% xlabel('Slice Encoding (ky)');
% ylabel('Phase Encoding (kz)');
% 
% % Setup the colorbar to correctly display the frame associations
% h_bar = colorbar;
% ylabel(h_bar, 'Time Frame (Region B) / Pass (Region A)');
% clim([min_Bj, max_Bj]); % Set color limits to the full data range
% 
% % --- 4. Animation Loop and Video Recording ---
% % Loop through every single point in the final, ordered sampling table.
% num_points_total = height(animation_table);
% for t = 1:num_points_total
%     % Get the linear index of the current point to sample
%     current_index = animation_table.Index(t);
% 
%     % Get the Bj value, which determines the point's color
%     current_bj = animation_table.Bj(t);
% 
%     % Place the point in our animation grid.
%     kspace_animation(current_index) = current_bj;
% 
%     % Update the 'CData' of the existing image plot for efficiency.
%     set(h_img, 'CData', kspace_animation);
% 
%     % *** UPDATE ALPHADATA IN THE LOOP ***
%     % As we fill in the kspace_animation matrix, we also need to update
%     % the AlphaData to make the newly filled pixels opaque.
%     set(h_img, 'AlphaData', ~isnan(kspace_animation));
% 
%     % Update the title to show the progress of the acquisition
%     title('k-space Sampling Animation');
% 
%     %Force MATLAB to draw the update.
%     drawnow; % Use 'drawnow' instead of 'limitrate' for smoother video capture
% 
%     %*** CAPTURE AND WRITE VIDEO FRAME ***
%     frame = getframe(gcf); % Capture the current figure window
%     writeVideo(v, frame);  % Write the captured frame to the video file
% 
%   if t>1 && t<num_points_total
%     if kspaceSamplingOrder_full.Bj(t) == (kspaceSamplingOrder_full.Bj(t-1)+1)
%       for ii = 1:100
%           writeVideo(v,frame)
%       end
%     end
%    end
% 
% end
% 
% title('k-space Sampling Animation - Finished');
% 
% %--- 5. Finalize Video ---
% close(v); % Close the video file to finalize and save it
% fprintf('Animation finished and video saved to: %s\n', videoFilename);