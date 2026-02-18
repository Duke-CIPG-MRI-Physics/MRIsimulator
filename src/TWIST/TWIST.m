function [TWIST_sampling_order] = TWIST(pA,pB,Matrix_Size_Acquired,FOV_acquired,R,PF_Factor)
%Roberto Carrascosa, Duke University, June 2025
% This function implements the TWIST sampling scheme for MRI
% according to:
%Song T, Laine AF, Chen Q, Rusinek H, Bokacheva L, Lim RP, Laub G,
%  Kroeker R, Lee VS. Optimal k-space sampling for dynamic
%  contrast-enhanced MRI with an application to MR renography. 
% Magn Reson Med. 2009 May;61(5):1242-8. doi: 10.1002/mrm.21901.
%  PMID: 19230014; PMCID: PMC2773550.

%This function outputs a table of points, in the order they should be
%sampled for TWIST. Note that TWIST operates only on 3D acquisitions.
%Therefore, actual use of this function will require further implemenetation 
%of frequency encoding, based on the user's needs.

%The inputs of this function are:
% pA --- defines the size of the central region that is sampled for
%   every measurement

% N --- the reciporal of pB, whre pB defines the fraction of the exterior
%  of k-space which is sampled every measurement. N must be an integer

% kspace_size --- vector of form: [#frequency,#phase (rows),#slice (columns)]



%The output of this function is a table of form:
% Linear Index, Bj, Row(phase), Column(slice)
%
% It is in the correct sampling order to achieve TWIST


%% --- Checking for validity of inputs
arguments
        % pA must be a single numeric value between 0.04 and 1, inclusive.
        pA (1,1) {mustBeNumeric, mustBeGreaterThanOrEqual(pA, 0), mustBeLessThanOrEqual(pA, 1)}

        % pB must be one of the allowed values (0,.1,.25,.33,.5)
        pB {mustBeMember(pB,[0,.1,.25,.33,.5])}

        % kspaceSize must be a 1x3 vector of positive, integer values.
        Matrix_Size_Acquired (1,3) {mustBeNumeric, mustBePositive, mustBeInteger}

        FOV_acquired (1,3) {mustBeNumeric, mustBePositive}

        %GRAPPA acceleration factors: [phase (rows), slice (columns)]
        R (1,2) {mustBeNumeric,mustBeInteger,mustBePositive} 
    
        %Partial Fourier acceleration factors: [phase (rows), slice (columns)]
        %Defaults to [1,1]
        PF_Factor (1,2) {mustBeNumeric, mustBePositive, mustBeLessThanOrEqual(PF_Factor,1), mustBeGreaterThan(PF_Factor,.5)}
end

Nb = round(1/pB);

%% --- Coordinate Grid Setup (Phase/Slice) --

kyi = 1:Matrix_Size_Acquired(3);  % slice (columns)
kzi = 1:Matrix_Size_Acquired(2);  % phase (rows)

centerPixel = floor(Matrix_Size_Acquired/2)+1;  % Use full kspaceSize as per your request
ky = kyi - centerPixel(3);  % col offset (slice)
kz = kzi - centerPixel(2);  % row offset (phase)
[kyM, kzM] = meshgrid(ky, kz);  % rows, cols

theta_matrix = cart2pol(kyM, kzM);

%Fix theta range, direction, and 0
theta_matrix(theta_matrix>0) = 2*pi-(theta_matrix(theta_matrix>0)); 
theta_matrix = abs(theta_matrix);


%% --- Get Region A---
[regionA,frequency_table] = getRegionA(Matrix_Size_Acquired,FOV_acquired,pA,PF_Factor,R);

%% --- Define Region B
regionB = ~regionA;

%% --- Radial and Angular Sorting ---
columnNames = ["Linear Index","Frequency","Theta","Region A?"];
unsortedData = table(frequency_table{:,1}, frequency_table{:,2}, theta_matrix(frequency_table{:,1}), regionA(frequency_table{:,1}),'VariableNames', columnNames);
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

bj_sequence = mod(0:(numRegionB_points - 1), Nb) + 1;

% Assign the generated 'bj_sequence' to the 'Bj' column for the
% identified Region B rows.
sortedData.Bj(regionB_rows) = bj_sequence';

sortedData_regionA = sortedData(sortedData.("Region A?"),:);
sortedData_regionB = sortedData(~sortedData.("Region A?"),:);



%% --- Creating Sampling Order

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

kspaceSamplingOrder_all_frames = [];
kspaceSamplingOrder_B = [];

for ii = 1:Nb
    B_current_Bj = sortedData_regionB(sortedData_regionB.Bj == ii,:);

    %Then in region B the trajectory Bj is sampled by first acquiring 
    % every odd point on the way outwards
    B_outward_current_Bj = B_current_Bj(1:2:end,:);
   
    %and then every even point on the way back in 
    B_inward_current_Bj = flipud(B_current_Bj(2:2:end,:));

    kspaceSamplingOrder_B_current_Bj = [B_outward_current_Bj ; B_inward_current_Bj];
    
    %This defines sampling order for only B region, as is used only in full
    %acquisition
    kspaceSamplingOrder_B = [kspaceSamplingOrder_B;kspaceSamplingOrder_B_current_Bj];


    %This creates the sampling order for the frames
    kspaceSamplingOrder_current_frame = [kspaceSamplingOrder_A;kspaceSamplingOrder_B_current_Bj];
    kspaceSamplingOrder_current_frame.("Frame") = ii * ones(height(kspaceSamplingOrder_current_frame),1);
    kspaceSamplingOrder_all_frames = [kspaceSamplingOrder_all_frames;kspaceSamplingOrder_current_frame];
    
end



%% --- Appending the initial full k-space acquisition

%assuming the overall trajectory is the same, but we only collect A at the
% very end

kspaceSamplingOrder_initial = [kspaceSamplingOrder_B;kspaceSamplingOrder_A];
kspaceSamplingOrder_initial.("Frame") = zeros(height(kspaceSamplingOrder_initial),1);


kspaceSamplingOrder_full = [kspaceSamplingOrder_initial;kspaceSamplingOrder_all_frames];


%% --- Cleaning up the output table

TWIST_sampling_order = kspaceSamplingOrder_full(:,{'Linear Index','Bj','Frame'});

[row,col] = ind2sub(Matrix_Size_Acquired(2:3),TWIST_sampling_order.("Linear Index"));
TWIST_sampling_order.("Row (phase)") = row;
TWIST_sampling_order.("Column (slice)") = col;

