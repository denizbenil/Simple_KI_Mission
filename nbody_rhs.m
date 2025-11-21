function [dxdt] = nbody_rhs(t, x, bodies, frame)
    % NBODY_RHS Evaluates the right-hand-side of a N-body propagator
    %
    % Inputs:
    %   t      - Current time [ephemeris seconds]
    %   x      - State vector [position; velocity] [km; km/s]
    %   bodies - Cell array of body structures
    %   frame  - Reference frame (e.g., 'ECLIPJ2000')
    %
    % Outputs:
    %   dxdt - Time derivative of state vector
    
    dxdt = zeros(6,1);
    dxdt(1:3) = x(4:6);  % velocity
    rr_ssb_obj = x(1:3);
    
    % Sum gravitational accelerations from all bodies
    for i = 1:length(bodies)
        rv_ssb_body = cspice_spkezr(bodies{i}.name, t, frame, 'NONE', 'Sun');
        rr_body_obj = rr_ssb_obj - rv_ssb_body(1:3);
        dist2 = dot(rr_body_obj, rr_body_obj);
        dist = sqrt(dist2);
        aa_grav = -bodies{i}.GM * rr_body_obj/(dist*dist2);
        dxdt(4:6) = dxdt(4:6) + aa_grav;
    end 
end