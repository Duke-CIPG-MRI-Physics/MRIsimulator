function [Contrast,outputArg2] = Ultrafast_Breast_Sim(Simulation_Parameters)
%Performs simulation of Ultrafast Breast MRI using analytically defined
%phantom
%   Detailed explanation goes here
arguments (Input)
    Simulation_Parameters struct
end

%% Input Deconstruction
scan_parameters = SimulationParameters.ScanParameters;
rBW_HzPerPix = SimulationParameters.ContrastParameters.rBW_HzPerPix;
TR = SimulationParameters.ContrastParameters.TR;
TE = SimulationParameters.ContrastParameters.TE;
pA = SimulationParameters.TWIST.pA;
pB = SimulationParameters.TWIST.pB;
Time_Measured = SimulationParameters.TWIST.Time_Measured;
R = SimulationParameters.ParallelImaging.GRAPPA_R;
PF_Factor = SimulationParameters.ParallelImaging.PF_Factor;

%% Setup
freq_phase_slice = [2 1 3]; % 1 = R/L, 2=A/P, 3 = S/I 
encodingFullStr = formatEncodingString(freq_phase_slice);

breastPhantomParams = createBreastPhantomParams();

% Derived contrast paramters
rBW_Hz = rBW_HzPerPix*matrix_size_acquired(1);  
dt_s = 1/rBW_Hz;   % dwell time between frequency-encode samples [s]

%% Configure acquisition ordering and timing
Sampling_Table = Ultrafast_Sampling(matrix_size_acquired,FOV_acquired,pA,Nb,Time_Measured,TR,R,PF_Factor);

k_idx_freq_pha_sli = [Sampling_Table.Frequency, Sampling_Table.("Row (phase)"), Sampling_Table.("Column (slice)")];

end