function EAxyz = quatToEAxyz(quat)

% Rot = Rx(EA(1)) * Ry(EA(2)) * Rz(EA(3))

R = quatToRot(quat);
EAxyz = rotToEAxyz(R);

end