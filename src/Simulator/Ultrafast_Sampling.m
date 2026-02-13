function [Sampling_Table,TWIST_Stats] = Ultrafast_Sampling(Matrix_Size_Acquired,FOV_acquired,pA,pB,N_Measurements,TR,R,PF_Factor)
%This function implements full sampling trajectory for Ultrafast MRI
%Imaging using TWIST,GRAPPA, and Partial Fourier
%   GRAPPA and PF are optional input. Due to MATLAB limitations, the way to
%   skip one or both inputs is to use a value of [1,1]

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

if all(R ~= 1)
    %GRAPPA
    Sampling_Table = GRAPPA_Undersample(Matrix_Size_Acquired,Sampling_Table,R);
end

if all(PF_Factor ~= 1)
    %Partial Fourier
    Sampling_Table = PF_Undersample(Matrix_Size_Acquired,Sampling_Table,PF_Factor);
end


%% --- Correcting for number of measurements 

Sampling_Table_A = Sampling_Table(Sampling_Table.Bj == 0, :);
Sampling_Table_B = Sampling_Table(Sampling_Table.Bj ~= 0, :);

Complete_Sampling_Table = Sampling_Table;
Complete_Sampling_Table.Bj(:) = 0;

frame_number = 1;
i_Bj = 1;
while frame_number <= N_Measurements 
    Frame_Sampling_Table = [Sampling_Table_A;Sampling_Table_B(Sampling_Table_B.Bj==i_Bj,:)];
    Complete_Sampling_Table = [Complete_Sampling_Table;Frame_Sampling_Table];
    
    frame_number = frame_number + 1;
    i_Bj = i_Bj + 1;

    if i_Bj > N
        i_Bj = 1;
    end
end

% % 2. Build the full periodic table in a single, vectorized operation
% if N_Measurements > 0 && ~isempty(Sampling_Table_B)
%     % Replicate the base table 'N_Measurements' times using repmat
%     Sampling_Table_B_Replicated = repmat(Sampling_Table_B, N_Measurements, 1);
% 
%     % Create a column vector of offsets to add (e.g., [0; N; 2*N; ...])
%     offsets = (0:N_Measurements-1)' * N;
% 
%     % Use repelem to expand the offsets to match the replicated table's size
%     bj_offsets = repelem(offsets, height(Sampling_Table_B), 1);
% 
%     % Add the offsets to the Bj column all at once
%     Sampling_Table_B_Replicated.Bj = Sampling_Table_B_Replicated.Bj + bj_offsets;
% 
%     % Remove rows that exceed the total number of measurements
%     Sampling_Table_B_Replicated(Sampling_Table_B_Replicated.Bj > Num_Measurements, :) = [];
% else
%     % Handle case where no repetitions are needed or Table B is empty
%     Sampling_Table_B_Replicated = Sampling_Table_B;
%     Sampling_Table_B_Replicated(Sampling_Table_B_Replicated.Bj > Num_Measurements, :) = [];
% end
% 
% % 3. Combine the static part with the new periodic part
% Sampling_Table = [Sampling_Table_A; Sampling_Table_B_Replicated];

%% --- Calculating Actual Time
Num_Measurements_Actual = Num_Measurements;
Temporal_Resolution_Actual = TR*sum(Sampling_Table.Bj ~= 0)/Num_Measurements_Actual;
Preparation_Scan_Time_Actual = TR*sum(Sampling_Table.Bj == 0);
Measurement_Time_Actual = Num_Measurements_Actual*Temporal_Resolution_Actual;

Timing_Actual = table(Temporal_Resolution_Actual,Num_Measurements_Actual,Preparation_Scan_Time_Actual,Measurement_Time_Actual);
Timing_Actual.Properties.VariableNames = {'Average Temporal Resolution (s)','# of Measurements','Preparation Scan Time (s)','Measurement Time (s)'};
Timing_Actual.Properties.RowNames = {'Actual Scan Timing'};

disp(Timing_Actual)

%% --- Incorporate Frequency Encoding
%This can take place at the very end as all we are doing is duplicating
%each entry by the number of Frequency encodes, this does not impact the
%calculation of scan time

N_freqs = Matrix_Size_Acquired(1); % We want to duplicate each row 3 times

% Step 1: Expand the table
% We generate a list of indices where each index is repeated N times
expandedIdx = repelem(1:height(Sampling_Table), N_freqs);

original_height = height(Sampling_Table);

Sampling_Table = Sampling_Table(expandedIdx, :);

% Step 2: Add the repetition counter
% We generate a vector [1; 2; 3] and stack it for every original row
Sampling_Table.Frequency = repmat((1:N_freqs)', original_height, 1);

%Correct linear index column for addded dimension
% Order: [Frequency, Phase, Slice, Frame]
sz_4D = [Matrix_Size_Acquired(1), Matrix_Size_Acquired(2), Matrix_Size_Acquired(3), max(Sampling_Table.Bj)+1];

Sampling_Table.("Linear Index") = sub2ind(sz_4D, ...
    Sampling_Table.Frequency, ...          % Dim 1
    Sampling_Table.("Row (phase)"), ...    % Dim 2
    Sampling_Table.("Column (slice)"), ... % Dim 3
    Sampling_Table.Bj + 1);                % Dim 4

end