function EA = rotToEulerAngles_new(R)

beta = atan2(-R(3,1),sqrt(R(1,1)^2+R(2,1)^2));   % [-90, 90]
alpha = atan2(R(2,1), R(1,1));
gamma = atan2(R(3,2), R(3,3));

EA = [alpha, beta, gamma]';

end