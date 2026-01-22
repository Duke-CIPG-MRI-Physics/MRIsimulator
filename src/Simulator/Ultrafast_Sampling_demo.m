clear;clc; close all;
%this demo illustrates that areas being sampled for Ultrafast MRI
%it incorporates most of the same logic as Ultrafast_Sampling.m


%% --- Setup
Complete_matrix_size = [224,224,224]; % [frequency, phase (rows), slice (columns)]

pA = .02; %Size of A Region for TWIST
N = 10; %Number of regions, Bj, to subdivide region B into

TR = 5.88E-3; %Sequence TR
Time_Measured = 90; %Time you wish to collect data for, excluding preparation scan time

GRAPPA_R = [3,2]; %GRAPPA acceleration factors: [phase (rows), slice (columns)]

Partial_Fourier_Factor = [6/8 , 6/8]; %Partial Fourier acceleration factors: [phase (rows), slice (columns)]

%% --- TWIST Sampling

Sampling_table = TWIST(pA,N,Complete_matrix_size);

Sampling_table_original = Sampling_table;

%% --- Timing Setup

%TWIST
Preparation_Scan_Time = TR*Complete_matrix_size(2)*Complete_matrix_size(3);
TWIST_Temporal_Resolution =  TR*(sum(Sampling_table.Bj ~= 0))/N;
Num_Measurements_TWIST = ceil(Time_Measured/TWIST_Temporal_Resolution);
Measurement_Time = TWIST_Temporal_Resolution*Num_Measurements_TWIST;

Timing = table(TWIST_Temporal_Resolution,Num_Measurements_TWIST,Preparation_Scan_Time,Measurement_Time);
Timing.Properties.VariableNames = {'Average Temporal Resolution (s)','# of Measurements','Preparation Scan Time (s)','Measurement Time (s)'};
Timing.Properties.RowNames = {'TWIST Only'};

%TWIST + GRAPPA
GRAPPA_acceleration_factor = 1/(GRAPPA_R(1)*GRAPPA_R(2));

Preparation_Scan_Time_GRAPPA = Preparation_Scan_Time*GRAPPA_acceleration_factor;
GRAPPA_Temporal_Resolution =  TWIST_Temporal_Resolution*GRAPPA_acceleration_factor;
Num_Measurements_GRAPPA = ceil(Time_Measured/GRAPPA_Temporal_Resolution);
Measurement_Time_GRAPPA = GRAPPA_Temporal_Resolution*Num_Measurements_GRAPPA;

Timing = [Timing ; {GRAPPA_Temporal_Resolution,Num_Measurements_GRAPPA,Preparation_Scan_Time_GRAPPA,Measurement_Time_GRAPPA}];
Timing.Properties.RowNames{end} = 'TWIST+GRAPPA';

%TWIST + Partial Fourier
PF_acceleration_factor = Partial_Fourier_Factor(1)*Partial_Fourier_Factor(2);

Preparation_Scan_Time_PF = Preparation_Scan_Time*PF_acceleration_factor;
PF_Temporal_Resolution =  TWIST_Temporal_Resolution*PF_acceleration_factor;
Num_Measurements_PF = ceil(Time_Measured/PF_Temporal_Resolution);
Measurement_Time_PF = PF_Temporal_Resolution*Num_Measurements_PF;

Timing = [Timing ; {PF_Temporal_Resolution,Num_Measurements_PF,Preparation_Scan_Time_PF,Measurement_Time_PF}];
Timing.Properties.RowNames{end} = 'TWIST+Partial Fourier';

%TWIST + GRAPPA + Partial Fourier
Preparation_Scan_Time_GRAPPA_PF = Preparation_Scan_Time_GRAPPA*PF_acceleration_factor;
GRAPPA_PF_Temporal_Resolution =  GRAPPA_Temporal_Resolution*PF_acceleration_factor;
Num_Measurements_GRAPPA_PF = ceil(Time_Measured/GRAPPA_PF_Temporal_Resolution);
Measurement_Time_GRAPPA_PF = GRAPPA_PF_Temporal_Resolution*Num_Measurements_GRAPPA_PF;

Timing = [Timing ; {GRAPPA_PF_Temporal_Resolution,Num_Measurements_GRAPPA_PF,Preparation_Scan_Time_GRAPPA_PF,Measurement_Time_GRAPPA_PF}];
Timing.Properties.RowNames{end} = 'TWIST+GRAPPA+Partial Fourier';

disp(Timing)
fprintf('\n***TIMES ARE APPROXIMATE***')
pause(1)

%% --- Asking what to use
GRAPPA_ask = input('\n\nUse GRAPPA? (y or n)\n','s')  ;
PF_ask = input('Use Partial Fourier? (y or n)\n','s')  ;

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

%% --- Integrating GRAPPA

if GRAPPA_ask == 'y'
fprintf('Removing GRAPPA Lines...\n')

%The strategy is to zero lines according to the R facor
Rows_to_keep = GRAPPA_R(1):GRAPPA_R(1):Complete_matrix_size(2);
Columns_to_keep = 1:GRAPPA_R(2):Complete_matrix_size(3);

Sampling_table_accelerated = [];

for ii = 1:height(Sampling_table)
    if any(Sampling_table{ii,"Row (phase)"} == Rows_to_keep) && any(Sampling_table{ii,"Column (slice)"} == Columns_to_keep)
       
        Sampling_table_accelerated = [Sampling_table_accelerated;Sampling_table(ii,:)];
 
    end

end

Sampling_table = Sampling_table_accelerated;

elseif GRAPPA_ask == 'n'
    fprintf('Skipping GRAPPA...\n')
end

%% --- Integrating Partial Fourier

if PF_ask == 'y'
fprintf('Removing PF Areas...\n')

%Deleting PF areas
Partial_Fourier_Number = [round(Complete_matrix_size(2)*Partial_Fourier_Factor(1)), round(Complete_matrix_size(3)*Partial_Fourier_Factor(2))];

Rows_to_keep = 1:Partial_Fourier_Number(1);
Columns_to_keep = (Complete_matrix_size(3)-Partial_Fourier_Number(2)):Complete_matrix_size(3);

Sampling_table_accelerated = [];

for ii = 1:height(Sampling_table)
    if any(Sampling_table{ii,"Row (phase)"} == Rows_to_keep) && any(Sampling_table{ii,"Column (slice)"} == Columns_to_keep)
       
        Sampling_table_accelerated = [Sampling_table_accelerated;Sampling_table(ii,:)];
 
    end

end

Sampling_table = Sampling_table_accelerated;

elseif PF_ask == 'n'
    fprintf('Skipping Partial Fourier...\n')
end

%% --- Correcting Number of Measurements 
Region_A_rows = (Sampling_table.Bj == 0);
Region_B_rows = (Sampling_table.Bj ~= 0);

Sampling_table_A = Sampling_table(Region_A_rows,:);
Sampling_table_B = Sampling_table(Region_B_rows,:);

Sampling_table_B_Corrected = Sampling_table_B;

for ii = 1:ceil(Num_Measurements/N)-1
    Sampling_table_B.Bj = Sampling_table_B.Bj+N;
    Sampling_table_B_Corrected = [Sampling_table_B_Corrected;Sampling_table_B];
end

Sampling_table_B_Corrected(Sampling_table_B_Corrected.Bj > Num_Measurements,:) = [];

Sampling_table = [Sampling_table_A;Sampling_table_B_Corrected];

%% --- Calculating Actual Time
Num_Measurements_Actual = Num_Measurements;
Temporal_Resolution_Actual = TR*sum(Sampling_table.Bj ~= 0)/Num_Measurements_Actual;
Preparation_Scan_Time_Actual = TR*sum(Sampling_table.Bj == 0);
Measurement_Time_Actual = Num_Measurements_Actual*Temporal_Resolution_Actual;

Timing_Actual = table(Temporal_Resolution_Actual,Num_Measurements_Actual,Preparation_Scan_Time_Actual,Measurement_Time_Actual);
Timing_Actual.Properties.VariableNames = {'Average Temporal Resolution (s)','# of Measurements','Preparation Scan Time (s)','Measurement Time (s)'};
Timing_Actual.Properties.RowNames = {'Actual Scan Timing'};

disp(Timing_Actual)

sampling_percentage = 100*(sum(Sampling_table.Bj ~= 0)/Num_Measurements_Actual)/(Complete_matrix_size(2)*Complete_matrix_size(3));

fprintf('\nEach frame is only sampling %.2f%% of k-space (on average)!!\n\n',sampling_percentage)

%% --- Displaying

filled_kSpace = zeros(Complete_matrix_size(2:3));
filled_kSpace(Sampling_table.("Linear Index")) = Sampling_table.Bj;

filled_kSpace_original = zeros(Complete_matrix_size(2:3));
filled_kSpace_original(Sampling_table_original.("Linear Index")) = Sampling_table_original.Bj;

figure
subplot(1,2,1)
imagesc(filled_kSpace)
colormap('copper')
title('Ultrafast Acquisition');

subplot(1,2,2)
imagesc(filled_kSpace_original)
colormap('copper')
title('Regular Acquisition');


filled_kSpace_frames = zeros([Complete_matrix_size(2:3),N]);

for ii = 1:N
    table_frame_rows = (Sampling_table.Bj == ii);
    frame_rows = Sampling_table.("Row (phase)")(table_frame_rows);
    frame_columns = Sampling_table.("Column (slice)")(table_frame_rows);

    for jj = 1:height(frame_rows)
        filled_kSpace_frames(frame_rows(jj),frame_columns(jj),ii) = 1;
    end
end

figure
imshow3D(filled_kSpace_frames)
title('Ultrafast Frames')

