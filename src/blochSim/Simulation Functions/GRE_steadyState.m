function [ Mss_postRF, Mss_TE, Mss_preRF ]  = GRE_steadyState(flip_angle_radians,T1_s,T2_s,TE_s,TR_s,dfreq_Hz)

% Calculate rotation matrix for RF pulse
Rflip = rotateAboutX(flip_angle_radians);

% Calculate a single matrix that simulates all the precession and decay
% between the RF pulse and the TE (where we measure our signal)
[Ate,Bte] = precessAndRelax(T1_s,T2_s,TE_s,dfreq_Hz);

% Calculate a single matrix that simulates all the precession and decay
% between the TE (where we measure our signal) and the TR (where our next
% RF pulse starts)
[Atr,Btr] = precessAndRelax(T1_s,T2_s,TR_s-TE_s,dfreq_Hz);

% Define:
% M0 - the magnetization immediately preceding the RF pulse
% M1 - the magnetization immediately following the RF pulse
% M2 - the magnetization TE
% M3 - the magnetization at TR
% 
% M1 = Rflip*M0
% M2 = Ate*M1 + Bte
% M3 = Atr*M2 + Btr
%
% We ultimately care about M2 (the echo time where we measure data)
% Define MSS_TE = M2
% MSS_TE = Ate*M1 + Bte
% Substituting in...
% MSS_TE = Ate*Rflip*M0+Bte 
% If we are in steady state, M3 = M0 
% MSS_TE = Ate*Rflip*M3+Bte
% Plugging in M3...
% MSS_TE = Ate*Rflip*(Atr*M2 + Btr) + Bte
% Remember M2 = MSS_TE
% MSS_TE = Ate*Rflip*(Atr*MSS_TE + Btr) + Bte
% Expanding...
% MSS_TE = Ate*Rflip*Atr*MSS_TE + Ate*Rflip*Btr + Bte
% Moving all the MSS_TE to the same side...
% MSS_TE - Ate*Rflip*Atr*MSS_TE = Ate*Rflip*Btr + Bte
% Factor out the common MSS_TE
% MSS_TE*(eye(3) - Ate*Rflip*Atr) = Ate*Rflip*Btr + Bte
% Invert the problem
% MSS_TE = ((eye(3) - Ate*Rflip*Atr)) \ (Ate*Rflip*Btr + Bte)

% MSS_TE = (eye(3) - Ate*Rflip*Atr) \ (Ate*Rflip*Btr + Bte);
Mss_TE = inv(eye(3) - Ate*Rflip*Atr) * (Ate*Rflip*Btr + Bte);

Mss_preRF = Atr*Mss_TE + Btr;

Mss_postRF = Rflip*Mss_preRF;
end