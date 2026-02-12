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
Time_Measured = 500; %sec

R = 1; %[2 3]
PF_Factor = 1; %[6/8 6/8]

%creating structure file
SimulationParameters.ScanParameters = scan_parameters;

SimulationParameters.ContrastParameters.rBW_HzPerPix = rBW_HzPerPix;
SimulationParameters.ContrastParameters.TR = TR;
SimulationParameters.ContrastParameters.TE = TE;

SimulationParameters.TWIST.pA = pA;
SimulationParameters.TWIST.pB = 1/Nb;
SimulationParameters.TWIST.Time_Measured = Time_Measured;

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







