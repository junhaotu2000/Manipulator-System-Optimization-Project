function model = make_model(inertias, masses)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
model = autoTree(3, 1, pi/2, 1); 
model.Xtree{end+1} = plux( eye(3), [0.1, 0, 0] );
end

