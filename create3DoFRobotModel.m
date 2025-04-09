% % Define each link's rotational inertia (using diagonal matrices)
% rotInertia = {diag([0.02, 0.02, 0.01]), diag([0.02, 0.02, 0.01]), diag([0.02, 0.02, 0.01])};
% 
% % Define link masses (vector with 3 elements)
% mass = [1.0, 1.0, 1.0];
% 
% % test the function
% model = create3DoFRobotModel(rotInertia, mass);
% disp(model.NB) 
% disp(model.parent) 
% for i=1:length(model.Xtree)
%     disp(size(model.Xtree{i}))
% end
% for i=1:length(model.I)
%     disp(size(model.I{i}))
% end
% %% Function Definitions
% 
function model = create3DoFRobotModel(rotInertia_arr, mass)
    % create3DoFRobotModel constructs a 3-DOF robot model using spatical v2.
    %
    % Inputs:
    %   rotInertia - Cell array of 3x3 rotational inertia matrices (one per link)
    %   mass       - 1x3 vector representing the mass of each link
    %
    % Output:
    %   model - Structure with fields:
    %             NB     : Number of joints (3)
    %             parent : Parent indices ([0, 1, 2])
    %             jtype  : Joint types (all 'Rz')
    %             Xtree  : Cell array of fixed transforms (6x6 matrices)
    %             I      : Cell array of spatial inertias (6x6 matrices)
    rotInertia = {diag(rotInertia_arr(1,:)),diag(rotInertia_arr(2,:)),diag(rotInertia_arr(3,:))};
    NB = 3;
    if length(mass) ~= NB || length(rotInertia) ~= NB
        error('mass and rotInertia must contain 3 elements each!');
    end
    
    L = [0.32, 0.32, 0.32, 0];                  % Link length (adjustable)
    
    twist_angles = [0, pi/2, pi/2];  % Twist angles for each joint
    
    model.NB = NB;
    model.parent = [0, 1, 2];   % Joint 1 is attached to the base; others connect sequentially
    model.jtype = cell(NB, 1);
    model.Xtree = cell(NB+1, 1); % NB fixed transforms plus one end-effector offset
    model.I = cell(NB, 1);
    
    model.appearance.base = {'tiles', [-1 1; -1 1; 0 0], 0.5};

    for i = 1:NB
        model.jtype{i} = 'Rz';  % All joints are revolute about the z-axis
        
        % Fixed transform: rotate about y-axis by twist angle then translate along x-axis by L
        R = roty(twist_angles(i));
        p = [0; 0; L(i)];
        model.Xtree{i} = plux(R, p);
        
        % Assume the center of mass is at the middle of the link (along x)
        com = [L(i)/2; 0; 0];
        model.I{i} = spatialInertia(mass(i), rotInertia{i}, com);
    end
    model.appearance.body{1} = {{'cyl', [0,0,.05;0,0,-.05],.05}};
    model.appearance.body{2} = {{'cyl', [0,0,.05;0,0,-.05],.05},{'cyl', [0,0,0;L(3),0,0],0.02}};
    model.appearance.body{3} = {{'cyl', [0,0,.05;0,0,-.05],.05},{'cyl', [0,0,0;0,L(4),0],0.02}};
    
    % End-effector offset: translation along x by 0.1 m
    model.Xtree{NB+1} = plux(eye(3), [0; L(4); 0]);

end

function I_sp = spatialInertia(m, J, p)
    % spatialInertia computes the 6x6 spatial inertia matrix.
    % The matrix is given by:
    %   I_sp = [ J + m*(S'*S),   m*S';
    %            m*S,           m*I3 ]
    % where S is the skew-symmetric matrix of p.
    S = skew(p);
    I_sp = [J + m*(S'*S), m*S'; 
            m*S,         m*eye(3)];
end

function S = skew(v)
    % skew returns the 3x3 skew-symmetric matrix of a 3x1 vector v.
    % For v = [v1; v2; v3], the output is:
    %    [  0   -v3   v2;
    %      v3    0   -v1;
    %     -v2   v1    0 ]
    S = [    0,  -v(3),  v(2);
          v(3),     0, -v(1);
         -v(2),  v(1),    0];
end

function R = roty(theta)
    % roty returns the 3x3 rotation matrix for a rotation about the y-axis.
    % Theta is in radians.
    R = [cos(theta),  0, sin(theta);
              0,      1,      0;
         -sin(theta), 0, cos(theta)];
end

function X = plux(R, p)
    % plux constructs a 6x6 spatial transform matrix from a rotation matrix and a translation vector.
    % Input:
    %   R - 3x3 rotation matrix
    %   p - 3x1 translation vector
    % Output:
    %   X - 6x6 spatial transform matrix
    X = [R, zeros(3); 
         skew(p)*R, R];
end