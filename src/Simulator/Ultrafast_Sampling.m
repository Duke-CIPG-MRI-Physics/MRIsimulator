function [Complete_Sampling_Table,TWIST_Stats] = Ultrafast_Sampling(Matrix_Size_Acquired,FOV_acquired,pA,pB,N_Measurements,TR,R,PF_Factor)
%This function implements full sampling trajectory for Ultrafast MRI
%Imaging using TWIST,GRAPPA, and Partial Fourier


%The inputs of this function are:

% Matrix_Size_Acquired --- vector of form: [#frequency, #phase (rows), #slice (columns)]

% FOV_acquired --- vector of form: [freq FOV, phase FOV, slice FOV] 
%   units do not matter

% pA --- defines the size of the central region that is sampled for
%   every measurement

% pB --- where pB defines the fraction of the exteriorvof k-space which is 
%   sampled every measurement. There is a limited set of allowed inputs

% N_Measurements --- the number of timepoints/frames after the initial full
% k-space acquisition to collect

% TR --- MRI repetion time, time required to collect a single line of
% k-space

% R --- GRAPPA acceleration factor of form [phase, slice]
%   can be skipped with value of 1

% PF_Factor --- Partial Fourier acceleration factor of form [phase fraction, slice fraction]
%   can be skipped with value of 1


%Compared to Timing Demo, this function does not output to command line,
%passes as output

%% --- Checking for validity of inputs
arguments (Input)

    %Complete_Matrix_Size is the final image size [frequency,phase,slice]
    Matrix_Size_Acquired (1,3) {mustBeNumeric, mustBePositive, mustBeInteger}

    FOV_acquired (1,3) {mustBeNumeric, mustBePositive}

    % pA defines size of TWIST 'A' region
    pA (1,1) {mustBeNumeric, mustBeGreaterThanOrEqual(pA, 0), mustBeLessThanOrEqual(pA, 1)}

    % pB must be one of the allowed values (0,.1,.25,.33,.5)
    pB {mustBeMember(pB,[0,.1,.25,.33,.5])}
    
    % N_Measurements is number of TWIST frames to acquire
    N_Measurements (1,1) {mustBePositive, mustBeInteger}

    %Desired TR for Sequence
    TR (1,1) {mustBeNumeric, mustBePositive}

    %GRAPPA acceleration factors: [phase (rows), slice (columns)]
    R (1,2) {mustBeNumeric,mustBeInteger,mustBePositive} 
    
    %Partial Fourier acceleration factors: [phase (rows), slice (columns)]
    %Defaults to [1,1]
    PF_Factor (1,2) {mustBeNumeric, mustBePositive, mustBeLessThanOrEqual(PF_Factor,1), mustBeGreaterThan(PF_Factor,.5)}

end

N = round(1/pB);



%% --- TWIST

Sampling_Table = TWIST(pA,pB,Matrix_Size_Acquired,FOV_acquired,R,PF_Factor);

%% --- Undersampling
% Strategy is to simply remove entries which correspond with the 2
% techniques

if any(R ~= 1)
    %GRAPPA
    Sampling_Table = GRAPPA_Undersample(Matrix_Size_Acquired,Sampling_Table,R);
end

if any(PF_Factor ~= 1)
    %Partial Fourier
    Sampling_Table = PF_Undersample(Matrix_Size_Acquired,Sampling_Table,PF_Factor);
end


%% --- Correcting for number of measurements 

Sampling_Table_initial = Sampling_Table(Sampling_Table.Frame == 0,:);

Sampling_Table_A = Sampling_Table(Sampling_Table.Frame == 0 & Sampling_Table.Bj == 0,:);
Sampling_Table_B = Sampling_Table(Sampling_Table.Frame == 0 & Sampling_Table.Bj ~= 0,:);



Complete_Sampling_Table = Sampling_Table_initial;

frame_number = 1;
i_Bj = 1;
while frame_number <= N_Measurements 
    Frame_Sampling_Table = [Sampling_Table_A;Sampling_Table_B(Sampling_Table_B.Bj==i_Bj,:)];
    Frame_Sampling_Table.Frame = frame_number * ones(height(Frame_Sampling_Table),1);

    Complete_Sampling_Table = [Complete_Sampling_Table;Frame_Sampling_Table];
    
    frame_number = frame_number + 1;
    i_Bj = i_Bj + 1;

    if i_Bj > N
        i_Bj = 1;
    end
end


%% --- Calculating Timing Stats
Preparation_Scan_Time = TR*sum(Sampling_Table.Frame == 0);
Measurement_Time = TR*sum(Complete_Sampling_Table.Frame ~= 0);
Temporal_Resolution = Measurement_Time/N_Measurements;

TWIST_Stats = table(Preparation_Scan_Time,Temporal_Resolution,Measurement_Time);
TWIST_Stats.Properties.VariableNames = {'Preparation Scan Time (s)','Average Temporal Resolution (s)','Measurement Time (s)'};

%% --- Incorporate Frequency Encoding
%This can take place at the very end as all we are doing is duplicating
%each entry by the number of Frequency encodes, this does not impact the
%calculation of scan time
% 
% N_freqs = Matrix_Size_Acquired(1); 
% 
% % Step 1: Expand the table
% % We generate a list of indices where each index is repeated N times
% expandedIdx = repelem(1:height(Complete_Sampling_Table), N_freqs);
% 
% original_height = height(Complete_Sampling_Table);
% 
% Complete_Sampling_Table = Complete_Sampling_Table(expandedIdx, :);
% 
% % Step 2: Add the repetition counter
% % We generate a vector [1; 2; 3] and stack it for every original row
% Complete_Sampling_Table.Frequency = repmat((1:N_freqs)', original_height, 1);
% 
% %Correct linear index column for addded dimension
% % Order: [Frequency, Phase, Slice, Frame]
% sz_4D = [Matrix_Size_Acquired(1), Matrix_Size_Acquired(2), Matrix_Size_Acquired(3), max(Complete_Sampling_Table.Bj)+1];
% 
% Complete_Sampling_Table.("Linear Index") = sub2ind(sz_4D, ...
%     Complete_Sampling_Table.Frequency, ...          % Dim 1
%     Complete_Sampling_Table.("Row (phase)"), ...    % Dim 2
%     Complete_Sampling_Table.("Column (slice)"), ... % Dim 3
%     Complete_Sampling_Table.Bj + 1);                % Dim 4

end