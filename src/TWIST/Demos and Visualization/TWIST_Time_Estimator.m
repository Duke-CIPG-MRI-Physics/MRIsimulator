function [A_time,B_time,frame_time] = TWIST_Time_Estimator(PE_plane_size,pA,N,TR)

A_points = (max(PE_plane_size)*pA)^2 * pi;
A_time = A_points*TR;

B_points = ((PE_plane_size(1)*PE_plane_size(2))-A_points)/N;
B_time = B_points*TR;

frame_time = A_time + B_time;

end