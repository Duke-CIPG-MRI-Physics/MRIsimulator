function Ry = rotateAboutY(yangle_radians)
% This function returns a rotation matrix to rotate by any angle yangle_radians about the Y-axis
Ry = [cos(yangle_radians) 0 sin(yangle_radians);
    0 1 0;
    -sin(yangle_radians) 0 cos(yangle_radians)];
end