quat = xout(1:4,:);
pc = xout(5:7,:);
spatialVel = xout(8:13,:);

angular = spatialVel(1:3,:);
linear = spatialVel(4:6,:);

plot(tout, pc)

