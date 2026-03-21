%Simple script for creating SimulatioParameters structure file for use with
%UltrafastBreastSim.m

clear;clc;

%parameters
load("breast_ultrafast_scan_parameters.mat")

rBW_HzPerPix = 570;
TR = (5.88E-3);  
TE = 2.63E-3;

pA = 0.05;
Nb = 10;
N_measurements = 20; %sec

R = 1; %[2 3]
PF_Factor = 1; %[6/8 6/8]



%creating structure file
SimulationParameters.ScanParameters = scan_parameters;

% Lesion parameters
SimulationParameters.LesionParameters.lesionArrivalDelay_s = 1;
SimulationParameters.LesionParameters.lesionWashinType = "instant";
SimulationParameters.LesionParameters.lesionWashoutType = "washout";
SimulationParameters.LesionParameters.lesionPeakEnhancement = 1.6;
SimulationParameters.LesionParameters.lesionBaselineDeltaIntensity = 0;


% MRI contrast parameters
SimulationParameters.MRIContrastParameters.rBW_HzPerPix = rBW_HzPerPix;
SimulationParameters.MRIContrastParameters.TR = TR;
SimulationParameters.MRIContrastParameters.TE = TE;

% Ultrafast parameters
SimulationParameters.TWIST.pA = pA;
SimulationParameters.TWIST.pB = 1/Nb;
SimulationParameters.TWIST.N_measurements = N_measurements;
SimulationParameters.TWIST.shareMode = "forward";
SimulationParameters.TWIST.shareMethod = "single_anchor";
SimulationParameters.TWIST.shareTieBreaker = "future";

SimulationParameters.ParallelImaging.GRAPPA_R = R;
SimulationParameters.ParallelImaging.PF_Factor = PF_Factor;


%Save Output
save_ask = input('Save output?: (y/n)','s');

if strcmpi(save_ask, 'y')
    fprintf('Input desired filename, file will be saved as <filename>.mat\n')
    filename = input(':','s');

    fprintf('Saving...\n')
    save(filename,'SimulationParameters')
    fprintf('File saved as %s.mat\n', filename);
else
    fprintf('Output not saved')
end







