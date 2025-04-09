clear; close all hidden;
rotInertia = [0.02, 0.02, 0.01;
              0.02, 0.02, 0.01;
              0.02, 0.02, 0.01];
mass = [1.0, 1.0, 1.0];
[time_grid, I_eff, torque_traj, angle_traj] = PlannerFunction(rotInertia, mass);
