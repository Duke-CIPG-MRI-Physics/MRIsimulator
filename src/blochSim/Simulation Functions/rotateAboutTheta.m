function Mtheta = rotateAboutTheta(alpha_radians,theta_radians)
% This function returns a rotation matrix for a rotation about the axis 
% defined by y=x*tan(phase_angle)
%
% We use a z-rotation to align the axis with the x-axis so we can do a
% simple x rotation, then rotate around z again to account for the theta
% angle
Mtheta = rotateAboutZ(theta_radians)*rotateAboutX(alpha_radians)*rotateAboutZ(-theta_radians);
end