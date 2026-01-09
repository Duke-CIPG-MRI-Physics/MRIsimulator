function [Mfp1,Mfp2]=precessAndRelax(T1_sec,T2_sec,t_sec,df_Hz)
%
%	Function simulates free precession and relaxation (T1/T2)
%	over a time interval T, given relaxation times T1 and T2
%	and off-resonance df.  Times are in s, off-resonance in Hz.
%
%   To get the result:
%   M = Mfp1*M0 + Mfp2
theta = 2*pi*df_Hz*t_sec;	% Resonant precession, radians.
E1 = exp(-t_sec/T1_sec);	% T1 
E2 = exp(-t_sec/T2_sec);    % T2

Mfp1 = [E2 0 0;0 E2 0;0 0 E1]*rotateAboutZ(theta);
Mfp2 = [0; 0; 1-E1];
end