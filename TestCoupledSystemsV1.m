clear; clc;

% Step 1: Run the mechanical subsystem script
run('MechOptimization');

% Now, after this script finishes, we have `results`, `x_opt`, `mechOut` in the workspace
% (Because the script put them in the base or local workspace.)

% Step 2: use mechOut to call 3D robot model
robotModel = create3DoFRobotModel(...
                    mechOut.D, ...
                    mechOut.t, ...
                    mechOut.L, ...
                    mechOut.linkMasses);

disp('--- 3-DOF Robot Model ---');
disp(robotModel);
