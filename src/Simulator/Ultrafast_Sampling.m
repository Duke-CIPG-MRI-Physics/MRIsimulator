function [Sampling_Table,PF_ask,GRAPPA_ask] = Ultrafast_Sampling(Matrix_Size_Acquired,FOV_acquired,pA,N,Time_Measured,TR,R,PF_Factor)
%This function implements full sampling trajectory for Ultrafast MRI
%Imaging using TWIST,GRAPPA, and Partial Fourier
%   GRAPPA and PF are optional input. Due to MATLAB limitations, the way to
%   skip one or both inputs is to use a value of [1,1]

%% --- Checking for validity of inputs
arguments (Input)

    %Complete_Matrix_Size is the final image size [frequency,phase,slice]
    Matrix_Size_Acquired (1,3) {mustBeNumeric, mustBePositive, mustBeInteger}

    FOV_acquired (1,3) {mustBeNumeric, mustBePositive}

    % pA defines size of TWIST 'A' region
    pA (1,1) {mustBeNumeric, mustBeGreaterThanOrEqual(pA, 0), mustBeLessThanOrEqual(pA, 1)}

    % N defines how many subregions, 'Bj', that 'B' region is divided into
    N (1,1) {mustBeNumeric, mustBePositive, mustBeInteger}
    
    %Time_measured must be a single, positive value in seconds
    Time_Measured (1,1) {mustBeNumeric, mustBePositive}

    %Desired TR for Sequence
    TR (1,1) {mustBeNumeric, mustBePositive}

    %GRAPPA acceleration factors: [phase (rows), slice (columns)]
    R (1,2) {mustBeNumeric,mustBeInteger,mustBePositive} 
    
    %Partial Fourier acceleration factors: [phase (rows), slice (columns)]
    %Defaults to [1,1]
    PF_Factor (1,2) {mustBeNumeric, mustBePositive, mustBeLessThanOrEqual(PF_Factor,1), mustBeGreaterThan(PF_Factor,.5)}

end

if N > ((Matrix_Size_Acquired(2)*Matrix_Size_Acquired(3)))/2
    error('N too large, must be smaller than (#phase*#slice)/2')
end

%% --- TWIST

Sampling_Table = TWIST(pA,N,Matrix_Size_Acquired,FOV_acquired,R,PF_Factor);

%% --- Timing Setup and Estimate

%TWIST
Preparation_Scan_Time = TR*Matrix_Size_Acquired(2)*Matrix_Size_Acquired(3);
TWIST_Temporal_Resolution =  TR*(sum(Sampling_Table.Bj ~= 0))/N;
Num_Measurements_TWIST = ceil(Time_Measured/TWIST_Temporal_Resolution);
Measurement_Time = TWIST_Temporal_Resolution*Num_Measurements_TWIST;

Timing = table(TWIST_Temporal_Resolution,Num_Measurements_TWIST,Preparation_Scan_Time,Measurement_Time);
Timing.Properties.VariableNames = {'Average Temporal Resolution (s)','# of Measurements','Preparation Scan Time (s)','Measurement Time (s)'};
Timing.Properties.RowNames = {'TWIST Only'};

if all(R ~= 1)
%TWIST + GRAPPA
GRAPPA_acceleration_factor = 1/(R(1)*R(2));

Preparation_Scan_Time_GRAPPA = Preparation_Scan_Time*GRAPPA_acceleration_factor;
GRAPPA_Temporal_Resolution =  TWIST_Temporal_Resolution*GRAPPA_acceleration_factor;
Num_Measurements_GRAPPA = ceil(Time_Measured/GRAPPA_Temporal_Resolution);
Measurement_Time_GRAPPA = GRAPPA_Temporal_Resolution*Num_Measurements_GRAPPA;

Timing = [Timing ; {GRAPPA_Temporal_Resolution,Num_Measurements_GRAPPA,Preparation_Scan_Time_GRAPPA,Measurement_Time_GRAPPA}];
Timing.Properties.RowNames{end} = 'TWIST+GRAPPA';
end

if all(PF_Factor ~= 1)
%TWIST + Partial Fourier
PF_acceleration_factor = PF_Factor(1)*PF_Factor(2);

Preparation_Scan_Time_PF = Preparation_Scan_Time*PF_acceleration_factor;
PF_Temporal_Resolution =  TWIST_Temporal_Resolution*PF_acceleration_factor;
Num_Measurements_PF = ceil(Time_Measured/PF_Temporal_Resolution);
Measurement_Time_PF = PF_Temporal_Resolution*Num_Measurements_PF;

Timing = [Timing ; {PF_Temporal_Resolution,Num_Measurements_PF,Preparation_Scan_Time_PF,Measurement_Time_PF}];
Timing.Properties.RowNames{end} = 'TWIST+Partial Fourier';
end

if all(R ~= 1) && all(PF_Factor ~= 1)
%TWIST + GRAPPA + Partial Fourier
Preparation_Scan_Time_GRAPPA_PF = Preparation_Scan_Time_GRAPPA*PF_acceleration_factor;
GRAPPA_PF_Temporal_Resolution =  GRAPPA_Temporal_Resolution*PF_acceleration_factor;
Num_Measurements_GRAPPA_PF = ceil(Time_Measured/GRAPPA_PF_Temporal_Resolution);
Measurement_Time_GRAPPA_PF = GRAPPA_PF_Temporal_Resolution*Num_Measurements_GRAPPA_PF;

Timing = [Timing ; {GRAPPA_PF_Temporal_Resolution,Num_Measurements_GRAPPA_PF,Preparation_Scan_Time_GRAPPA_PF,Measurement_Time_GRAPPA_PF}];
Timing.Properties.RowNames{end} = 'TWIST+GRAPPA+Partial Fourier';
end
fprintf('TIMING ESTIMATES:\n')
disp(Timing)
fprintf('\n')


%% --- Asking what to use
if any(R ~= 1) && all(PF_Factor == 1)
PF_ask = 'n';

GRAPPA_ask = input('\n\nUse GRAPPA? (y or n)\n','s');

    if GRAPPA_ask == 'n' 
        Num_Measurements = Num_Measurements_TWIST;
    elseif GRAPPA_ask =='y'
        Num_Measurements = Num_Measurements_GRAPPA;

    else
        error('At least one input not supported, please respond y or n')
    end



elseif any(PF_Factor ~= 1) && all(R == 1)
GRAPPA_ask = 'n';

PF_ask = input('Use Partial Fourier? (y or n)\n','s');

    if PF_ask =='n'
        Num_Measurements = Num_Measurements_TWIST;
    elseif PF_ask =='y'
        Num_Measurements = Num_Measurements_PF;
    else
    error('At least one input not supported, please respond y or n')
    end

elseif any(R ~= 1) && any(PF_Factor ~= 1)
    
GRAPPA_ask = input('\n\nUse GRAPPA? (y or n)\n','s');
PF_ask = input('Use Partial Fourier? (y or n)\n','s');

    if GRAPPA_ask == 'n' && PF_ask =='n'
        Num_Measurements = Num_Measurements_TWIST;
    elseif GRAPPA_ask =='y' && PF_ask == 'n'
        Num_Measurements = Num_Measurements_GRAPPA;
    elseif GRAPPA_ask == 'n' && PF_ask == 'y'
        Num_Measurements = Num_Measurements_PF;
    elseif GRAPPA_ask == 'y' && PF_ask =='y'
        Num_Measurements = Num_Measurements_GRAPPA_PF;
    else
    error('At least one input not supported, please respond y or n')
    end
else
    GRAPPA_ask = 'n';
    PF_ask = 'n';
    Num_Measurements = Num_Measurements_TWIST;
end
%% --- Integrating GRAPPA
if GRAPPA_ask == 'y'

fprintf('Removing GRAPPA Lines...\n')
Sampling_Table = GRAPPA_Undersample(Matrix_Size_Acquired,Sampling_Table,R);

elseif PF_ask == 'n'
    fprintf('Skipping GRAPPA...\n')
    pause(.2)
end

%% --- Integrating Partial Fourier
if PF_ask == 'y'
   
fprintf('Removing PF Areas...\n')
Sampling_Table = PF_Undersample(Matrix_Size_Acquired,Sampling_Table,PF_Factor);

elseif PF_ask == 'n'
    fprintf('Skipping Partial Fourier...\n')
    pause(.2)
end

%% --- Correcting Number of Measurements 

% 1. Split the table into a static part (A) and a part to be repeated (B)
Region_A_rows = (Sampling_Table.Bj == 0);
Sampling_Table_A = Sampling_Table(Region_A_rows, :);
Sampling_Table_B = Sampling_Table(~Region_A_rows, :);

% 2. Calculate the total number of repetitions needed
num_reps = ceil(Num_Measurements / N);

% 3. Build the full periodic table in a single, vectorized operation
if num_reps > 0 && ~isempty(Sampling_Table_B)
    % Replicate the base table 'num_reps' times using repmat
    Sampling_Table_B_Corrected = repmat(Sampling_Table_B, num_reps, 1);
    
    % Create a column vector of offsets to add (e.g., [0; N; 2*N; ...])
    offsets = (0:num_reps-1)' * N;
    
    % Use repelem to expand the offsets to match the replicated table's size
    bj_offsets = repelem(offsets, height(Sampling_Table_B), 1);
    
    % Add the offsets to the Bj column all at once
    Sampling_Table_B_Corrected.Bj = Sampling_Table_B_Corrected.Bj + bj_offsets;
    
    % Remove rows that exceed the total number of measurements
    Sampling_Table_B_Corrected(Sampling_Table_B_Corrected.Bj > Num_Measurements, :) = [];
else
    % Handle case where no repetitions are needed or Table B is empty
    Sampling_Table_B_Corrected = Sampling_Table_B;
    Sampling_Table_B_Corrected(Sampling_Table_B_Corrected.Bj > Num_Measurements, :) = [];
end

% 4. Combine the static part with the new periodic part
Sampling_Table = [Sampling_Table_A; Sampling_Table_B_Corrected];

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
% Order: [Frequency, Phase, Slice, Time]
sz_4D = [Matrix_Size_Acquired(1), Matrix_Size_Acquired(2), Matrix_Size_Acquired(3), max(Sampling_Table.Bj)+1];

Sampling_Table.("Linear Index") = sub2ind(sz_4D, ...
    Sampling_Table.Frequency, ...          % Dim 1
    Sampling_Table.("Row (phase)"), ...    % Dim 2
    Sampling_Table.("Column (slice)"), ... % Dim 3
    Sampling_Table.Bj + 1);                % Dim 4

end
