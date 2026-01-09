function [MTrans, MLong] = calculateMxyAndMz(M)
    % This function takes inputs of a 3xn matrix, where each column represents a 3x1 vector in 3D. 
    % The function returns two 1xn matrices that correspond to the transverse magnetization (Mxy) and the longitudinal magnetization (Mz) respectively.
    MTrans = abs(M(1,:)+M(2,:)*1i);
    MLong = M(3,:);
end