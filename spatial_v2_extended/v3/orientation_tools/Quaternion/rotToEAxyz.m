function EAxyz = rotToEAxyz(R)

if isa(R(1,1), 'casadi.MX')
    EAxyz = casadi.MX(3,1);
else
    EAxyz = zeros(3,1);
end
EAxyz(2) = atan2(R(1,3),sqrt(R(1,1)^2+R(1,2)^2));   % [-90, 90]
EAxyz(3) = atan2(-R(1,2), R(1,1));
EAxyz(1) = atan2(-R(2,3), R(3,3));

end