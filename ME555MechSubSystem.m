%% runGA.m
clear; clc;
%rng(1234);  % Fix the random seed for reproducibility

%% 1. Variable Definition
nLinks = 3;  
% Design variable vector x = [D1,...,D3, t1,...,t3, L1,...,L3, m]
% D: Link outer diameter (mm), t: Wall thickness (mm), L: Link length (mm), m: Material selection (integer: 1~4)
nCont = nLinks*3;      % Number of continuous variables (D, t, L)
nVar = nCont + 1;      % Total number of variables
popSize = 300;         % Population size
nCrit = 50;            % Number of critical population individuals
numRuns = 5;           % Number of optimization runs

% Define the lower and upper bounds for design variables
lb = [ones(1, nLinks)*20,  ones(1, nLinks)*5,  ones(1, nLinks)*50, 1];
ub = [ones(1, nLinks)*200, ones(1, nLinks)*10, ones(1, nLinks)*800, 4];

% Specify the index of integer variables (the last variable represents material selection)
IntCon = nVar;

% Construct an initial population matrix
initialPop = zeros(popSize, nVar);

% Initialize critical population individuals with midpoint values
for i = 1:nCrit
    initialPop(i,:) = (lb + ub) / 2;
end

% Generate remaining individuals randomly
for i = (nCrit+1):popSize
    initialPop(i,:) = lb + rand(1, nVar) .* (ub - lb);
end

%% 2. Set GA Optimization Options
optionsGA = optimoptions('ga',...
    'Display','iter',...
    'PopulationSize',popSize,...
    'MaxGenerations',600,...
    'PlotFcn', {@gaplotbestf, @gaplotstopping},...
    'InitialPopulationMatrix', initialPop);

%% 3. Execute GA Optimization Multiple Times and Collect Data
results = cell(numRuns, 5);  % Table storage

for i = 1:numRuns
    fprintf('Starting GA optimization - Run %d...\n', i);
    
    % Execute GA
    [x_opt, fval_opt] = ga(@massObjective, nVar, [], [], [], [], lb, ub, @nonlcon, IntCon, optionsGA);
    
    % Compute constraint values at optimal point
    [c_active, ~] = nonlcon(x_opt);
    
    % Identify active constraints (close to zero)
    active_constraints = abs(c_active) < 1e-5;
    
    % Store results in table format
    results{i,1} = x_opt;               % Optimal design variables
    results{i,2} = fval_opt;            % Optimal objective value (mass)
    results{i,3} = initialPop(1,:);     % Representative initial point
    % results{i,4} = c_active;            % Constraint values
    % results{i,5} = active_constraints;  % Active constraints (boolean)
    
    fprintf('Run %d completed.\n', i);
end

%% 4. Convert Results to Table and Display
T = table(...
    results(:,1), results(:,2), results(:,3), ...
    'VariableNames', {'x_opt', 'f_opt', 'InitialPoint'});

% Display table
disp(T);


function f = massObjective(x)

nLinks = 3;
D = x(1:nLinks);                % Link outer diameter (mm)
t = x(nLinks+1:2*nLinks);       % Wall thickness (mm)
L = x(2*nLinks+1:3*nLinks);     % Link length (mm)
m = round(x(end));              % Material selection (rounded to integer)

% Material properties: Density (kg/m^3)
switch m
    case 1  % Aluminum alloy
        rho = 2700;
    case 2  % Carbon fiber
        rho = 1600;
    case 3  % Titanium alloy
        rho = 4500;
    case 4  % High-strength steel
        rho = 7800;
    otherwise
        rho = 2700;
end

% Compute total mass: Convert dimensions from mm to m
f = 0;
for i = 1:nLinks
    D_m = D(i) * 0.001;    % Convert to meters
    t_m = t(i) * 0.001;
    L_m = L(i) * 0.001;
    
    % Compute cross-sectional area of hollow cylinder (m^2)
    A = (pi/4)*((D_m)^2 - (D_m - 2*t_m)^2);
    V = A * L_m;      % Link volume (m^3)
    mass_i = rho * V; % Link mass (kg)
    f = f + mass_i;
end

f = 10*f + 20;  % Regulating and assuming additional 20kg for joints and other components

end

function [c, ceq] = nonlcon(x)

nLinks = 3;
D = x(1:nLinks);              % Link outer diameter (mm)
t = x(nLinks+1:2*nLinks);     % Wall thickness (mm)
L = x(2*nLinks+1:3*nLinks);   % Link length (mm)
m = round(x(end));            % Material selection

% Safety factor
S = 4;      
% Assumed maximum bending moment per link (Nm)
M_max = 500; 

% Material yield strength (MPa) and stiffness (GPa)
switch m
    case 1  % Aluminum alloy
        sigma_yield = 250e6;
        E = 70e9;
    case 2  % Carbon fiber
        sigma_yield = 500e6;
        E = 60e9;
    case 3  % Titanium alloy
        sigma_yield = 900e6;
        E = 110e9;
    case 4  % High-strength steel
        sigma_yield = 600e6;
        E = 210e9;
    otherwise
        sigma_yield = 250e6;
        E = 70e9;
end

c = [];  % Inequality constraints: c(x) <= 0

for i = 1:nLinks
    % Retrieve parameters for link i
    D_i = D(i);  % mm
    t_i = t(i);  % mm
    L_i = L(i);
    
    % Convert units to meters
    D_m = D_i * 0.001;
    t_m = t_i * 0.001;
    L_m = L_i * 0.001;
    
    % Moment of inertia I (m^4) for a hollow circular section
    I = (pi/64) * ((D_m)^4 - (D_m - 2*t_m)^4);
    
    % Compute maximum bending stress (Pa)
    sigma = (M_max * (D_m/2)) / I;
    
    % Strength constraint: sigma <= sigma_yield/S
    c(end+1) = sigma - (sigma_yield / S);
    
    % Manufacturing constraint: Inner diameter (D - 2*t) > 0, rewritten as 2*t - D <= 0
    c(end+1) = 2*t_i - D_i;

    % Deflection δ = (F_max * L^3)/(3*E*I)
    delta = (M_max * pi*((D_m)^2 - (D_m - 2*t_m)^2)/4 * L_m^3) / (3 * E * I);
    
    % Stiffness constraint: δ <= 0.005*L is required, that is, δ - 0.005*L_m <= 0
    c(end+1) = delta - 0.001 * L_m;


end

% Workspace constraint: Total link length must be >= 1000 mm
c(end+1) = 1000 - sum(L);

ceq = [];  % No equality constraints

end
