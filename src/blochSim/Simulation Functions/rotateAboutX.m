function Rx = rotateAboutX(xangle_radians)
% This function returns a rotation matrix to rotate by any angle xangle_radians about the X-axis
Rx = [1 0 0; 
    0 cos(xangle_radians) -sin(xangle_radians);
    0 sin(xangle_radians) cos(xangle_radians)];
end