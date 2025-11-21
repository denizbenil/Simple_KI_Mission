function [bodies] = nbody_init(labels)
    % NBODY_INIT Initialize planetary data for n-body propagation
    %
    % Inputs:
    %   labels - Cell array of body names (SPICE convention)
    %
    % Outputs:
    %   bodies - Cell array of structures with name and GM
    
    bodies = cell(size(labels));
    for i = 1:length(labels)
        bodies{i}.name = labels(i);
        bodies{i}.GM = cspice_bodvrd(labels{i},'GM',1);
    end
end
