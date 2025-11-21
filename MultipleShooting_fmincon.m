clc; clear all; close all; cspice_kclear();
format long;

% =========================================================================
% PLANETARY DEFENSE MISSION OPTIMIZATION FOR ASTEROID APOPHIS
% With Parametric Time Variables (alpha0, alpha1, TOF)
% =========================================================================

%% SPICE KERNEL LOADING
cspice_furnsh('kernels\naif0012.tls'); % (LSK)
cspice_furnsh('kernels\de432s.bsp');   % (SPK)
cspice_furnsh('kernels\20099942_Apophis.bsp');   % (SPK)
cspice_furnsh('kernels\gm_de432.tpc'); % (PCK)
cspice_furnsh('kernels\pck00010.tpc'); % (PCK)

fprintf('Number of LSK  kernels: %d\n', cspice_ktotal('lsk'));
fprintf('Number of SPK  kernels: %d\n', cspice_ktotal('spk'));
fprintf('Number of PCK  kernels: %d\n', cspice_ktotal('pck'));
fprintf('Number of CK   kernels: %d\n', cspice_ktotal('ck'));
fprintf('Number of TEXT kernels: %d\n', cspice_ktotal('TEXT'));
fprintf('\nTOTAL kernels number: %d\n', cspice_ktotal('ALL'));

%% INITIAL PARAMETERS
fprintf('\n=== INITIAL PARAMETERS ===\n');
Mast = 4e10; % Asteroid Apophis mass [kg]
betaparameter = 3; % Beta parameter (momentum enhancement factor)
Isp = 300; % Specific impulse [s]
g0 = 9.80665; % Standard gravity [m/s^2]

fprintf('Asteroid mass (Mast): %g kg\n', Mast);
fprintf('Beta parameter: %g\n', betaparameter);
fprintf('Specific Impulse (Isp): %g s\n', Isp);
fprintf('Standard gravity (g0): %g m/s^2\n', g0);
fprintf('NOTE: Momentum coefficient will be calculated dynamically based on final spacecraft mass\n');
fprintf('=========================================\n\n');

%% NATURAL TRAJECTORY ANALYSIS
% Define list of celestial bodies
labels = {'Sun';
          'Mercury';
          'Venus';
          'Earth';
          'Moon';
          'Mars Barycenter';
          'Jupiter Barycenter';
          'Saturn Barycenter';
          'Uranus Barycenter';
          'Neptune Barycenter';
          'Pluto Barycenter'};
bodies = nbody_init(labels);

% Select integration frame
center = 'SSB';
frame = 'ECLIPJ2000'; % ecliptic plane

% Asteroid natural trajectory
asteroid_label = '20099942';
ref_epoch_str = '2028-Sep-01 00:00:00.0000 TDB';
et0 = cspice_str2et(ref_epoch_str);
x0 = cspice_spkezr('20099942',et0,frame,'NONE',center);
final_epoch_str = '2029-May-1 00:00:00.0000 TDB';
etf = cspice_str2et(final_epoch_str);

tt = linspace(et0,etf,10000); 
xx = cspice_spkezr('20099942',tt,frame,'NONE',center);    
epoch = datetime('2000-01-01 12:00:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
dateTimes = epoch + seconds(tt);

% Get planetary positions
others = {'Earth', 'Moon', 'Sun'};
planet_colors = {'#77AC30','#A2142F','#EDB120'};
for i = 1:length(others)    
    rr_planet(:,:,i) = cspice_spkezr(others{i},tt,frame,'NONE',center);
end

% Calculate distances
EA_vec = rr_planet(1:3,:,1) - xx(1:3,:);
dist_EA = (EA_vec(1,:).^2 + EA_vec(2,:).^2 + EA_vec(3,:).^2).^0.5;
au_ea = cspice_convrt(dist_EA,'km','au');

MA_vec = rr_planet(1:3,:,2) - xx(1:3,:);
dist_MA = (MA_vec(1,:).^2 + MA_vec(2,:).^2 + MA_vec(3,:).^2).^0.5;
au_ma = cspice_convrt(dist_MA,'km','au');

% Find natural closest approach
[min_dist,ind]=min(dist_EA);
t_TCA_1=cspice_et2utc(tt(ind),'C',4);
t12h = linspace(tt(ind)-6*60*60,tt(ind)+6*60*60,101);

fprintf('\n=== APOPHIS NATURAL TRAJECTORY INFO ===\n');
fprintf('Natural TCA: %s\n', t_TCA_1);
fprintf('Min Distance: %.2f km (%.4f AU)\n', min_dist, cspice_convrt(min_dist,'km','au'));
fprintf('========================================\n\n');

%% OPTIMIZATION SETUP WITH PARAMETRIC TIME VARIABLES

d2s = 24*60*60;  % days to seconds
s2d = 1/d2s;     % seconds to days

% Define the time boundaries for APOPHIS mission
LWO     = cspice_str2et('2024-Jun-01 00:00:00.0000 TDB');
LWC     = cspice_str2et('2024-Dec-01 00:00:00.0000 TDB');
IMPS    = cspice_str2et('2028-Dec-01 00:00:00.0000 TDB');
IMPF    = cspice_str2et('2029-Apr-01 00:00:00.0000 TDB'); 

% Change reference for optimization
center = 'Sun';
frame = 'ECLIPJ2000';
labels = {'Sun'};
bodies = nbody_init(labels);

% Time bounds structure to pass to functions
time_bounds = struct('LWO', LWO*s2d, 'LWC', LWC*s2d, ...
                     'IMPS', IMPS*s2d, 'IMPF', IMPF*s2d);

% Initial guess for parametric time variables
alpha0_init = 0.3;     % Fraction of available launch window [0, 0.9]
alpha1_init = 0.5;     % Fraction of TOF for DSM timing [0, 0.9]
TOF_init = 1000;        % Total time of flight [days]

% Calculate derived times from parametric variables
Delta_t0 = (IMPF - LWO)*s2d - TOF_init;
t_launch_init = (LWO*s2d + Delta_t0 * alpha0_init);
t_IMP_init = t_launch_init + TOF_init;
t_DSM_init = t_launch_init + alpha1_init * TOF_init;

% Get initial guess for x2 and x3
rv_apophis_dsm = cspice_spkezr('20099942', t_DSM_init*d2s, frame, 'NONE', center);
rv_apophis_impact = cspice_spkezr('20099942', t_IMP_init*d2s, frame, 'NONE', center);

% Initial parametric launch variables
Delta_v0_init = 3.0;  % km/s - hyperbolic excess velocity magnitude
alpha_Delta_v0_init = pi/4;  % radians - in-plane angle
delta_Delta_v0_init = 0.1;   % radians - out-of-plane angle
msc0_init = 1000;  % kg - initial spacecraft mass

% X0 = [alpha0; alpha1; TOF; Delta_v0; alpha_Delta_v0; delta_Delta_v0; msc0; x2(pos,vel); x3(pos,vel)]
X0 = [alpha0_init; alpha1_init; TOF_init; ...                               % time params (3)
      Delta_v0_init; alpha_Delta_v0_init; delta_Delta_v0_init; msc0_init;  % launch+mass (4)
      rv_apophis_dsm(1:3); rv_apophis_dsm(4:6) + [5; 5; 0.001];                    % x2 (6)
      rv_apophis_impact(1:3); rv_apophis_impact(4:6)];                             % x3 (6)
% Total: 19 variables

% Constraint matrices and bounds
A = [];
b = [];
Aeq = [];
beq = [];

Rl = 1.4*10^9; % Saturn distance from sun

% Lower bounds: [alpha0, alpha1, TOF, Delta_v0, alpha_Delta_v0, delta_Delta_v0, msc0, x2(6), x3(6)]
lb = [0, 0.05, 1000, ...                                    % time params
      1.5, 0, 0, 100, ...                               % launch+mass
      -Rl, -Rl, -Rl, -100, -100, -100, ...             % x2
      -Rl, -Rl, -Rl, -100, -100, -100];                % x3

% Upper bounds
ub = [0.9, 0.95, 2000, ...                               % time params
      5.5, 2*pi, pi/2, 3000, ...                        % launch+mass
      Rl, Rl, Rl, 100, 100, 100, ...                   % x2
      Rl, Rl, Rl, 100, 100, 100];                      % x3

% Nonlinear constraint function & objective function definition
nonlcon = @(x)nonlinconstMultipleShoot_parametric_alpha(...
    x(1), x(2), x(3), ...           % alpha0, alpha1, TOF
    x(4), x(5), x(6), x(7), ...     % Delta_v0, alpha, delta, msc0
    x(8:13), x(14:19), ...          % x2, x3
    betaparameter, Mast, Isp, g0, time_bounds);

fun_obj = @(x)objFunMultipleShoot_parametric_alpha(...
    x(1), x(2), x(3), ...           % alpha0, alpha1, TOF
    x(4), x(5), x(6), x(7), ...     % Delta_v0, alpha, delta, msc0
    x(8:13), x(14:19), ...          % x2, x3
    betaparameter, Mast, Isp, g0, time_bounds);

% Optimization options
options = optimoptions('fmincon','SpecifyObjectiveGradient',false,...
    'SpecifyConstraintGradient',false,'Display','iter',...
    'Algorithm','active-set',MaxFunctionEvaluations=20000,...
    MaxIterations=10000,CheckGradients=false);

fprintf('\n=== STARTING OPTIMIZATION WITH PARAMETRIC TIME VARIABLES ===\n');
fprintf('Decision Variables:\n');
fprintf('  - alpha0: Launch time fraction [0, 0.9]\n');
fprintf('  - alpha1: DSM time fraction [0, 0.9]\n');
fprintf('  - TOF: Total time of flight [200, 1000] days\n');
fprintf('  - Delta_v0, angles, mass: Launch parameters\n');
fprintf('  - x2, x3: Trajectory states\n');
fprintf('Total: 19 optimization variables\n\n');

tic
[X_sol, fval, exitflag, output] = fmincon(fun_obj, X0, A, b, Aeq, beq, lb, ub, nonlcon, options);
toc

fprintf('\n=== OPTIMIZATION COMPLETED ===\n');
fprintf('Exit Flag: %d\n', exitflag);
fprintf('Function Evaluations: %d\n', output.funcCount);
fprintf('Iterations: %d\n', output.iterations);
fprintf('Final Objective Value: %.6f km (negative distance)\n', fval);
fprintf('==============================\n\n');

%% POST-PROCESSING
postprocess_parametric_alpha(X_sol, betaparameter, Mast, Isp, g0, time_bounds);

%% ========================================================================
%  FUNCTION DEFINITIONS
%  ========================================================================

function [t0, t_DSM, t_IMP] = parametric_to_times(alpha0, alpha1, TOF, time_bounds)
    % Convert parametric time variables to absolute ephemeris times
    %
    % Inputs:
    %   alpha0      - Launch time fraction [0, 0.9]
    %   alpha1      - DSM time fraction [0, 0.9]
    %   TOF         - Total time of flight [days]
    %   time_bounds - Structure with LWO, LWC, IMPS, IMPF [days]
    %
    % Outputs:
    %   t0      - Launch epoch [days]
    %   t_DSM   - DSM epoch [days]
    %   t_IMP   - Impact epoch [days]
    
    Delta_t0 = time_bounds.IMPF - time_bounds.LWO - TOF;
    t0 = time_bounds.LWO + Delta_t0 * alpha0;
    t_DSM = t0 + alpha1 * TOF;
    t_IMP = t0 + TOF;
end

function [r_departure, v_departure] = parametric_to_cartesian(Delta_v0, alpha_Delta_v0, delta_Delta_v0, t_launch, center, frame)
    % Convert parametric launch variables to Cartesian departure state
    %
    % Inputs:
    %   Delta_v0        - Magnitude of hyperbolic excess velocity [km/s]
    %   alpha_Delta_v0  - In-plane angle from tangent direction [rad]
    %   delta_Delta_v0  - Out-of-plane angle from orbital plane [rad]
    %   t_launch        - Launch epoch [ephemeris time in seconds]
    %   center          - Reference center (e.g., 'Sun')
    %   frame           - Reference frame (e.g., 'ECLIPJ2000')
    %
    % Outputs:
    %   r_departure - Departure position in inertial frame [km]
    %   v_departure - Departure velocity in inertial frame [km/s]
    
    % Step 1: Get Earth's state at launch
    rv_earth = cspice_spkezr('Earth', t_launch, frame, 'NONE', center);
    r_earth = rv_earth(1:3);
    v_earth = rv_earth(4:6);
    
    % Step 2: Construct v_infinity in TNH frame
    % T: Tangent (along velocity direction)
    % N: Normal (perpendicular to velocity in orbital plane)
    % H: Angular momentum direction (out of orbital plane)
    v_TNH = [Delta_v0 * cos(alpha_Delta_v0) * cos(delta_Delta_v0);
             Delta_v0 * sin(alpha_Delta_v0) * cos(delta_Delta_v0);
             Delta_v0 * sin(delta_Delta_v0)];
    
    % Step 3: Build TNH to inertial frame transformation matrix
    % T: tangent (along velocity)
    T_hat = v_earth / norm(v_earth);
    
    % H: angular momentum direction
    H_vec = cross(r_earth, v_earth);
    H_hat = H_vec / norm(H_vec);
    
    % N: completes right-handed system
    N_hat = cross(H_hat, T_hat);
    
    % Rotation matrix from TNH to inertial frame
    R_TNH_to_inertial = [T_hat, N_hat, H_hat];
    
    % Step 4: Transform v_infinity to inertial frame and add Earth velocity
    v_inf_inertial = R_TNH_to_inertial * v_TNH;
    v_departure = v_earth + v_inf_inertial;
    r_departure = r_earth;
end

function [value, isterminal, direction] = closest(t, x, center, frame)
    % Event function to detect closest approach to Earth
    x_E = cspice_spkezr('Earth', t, frame, 'NONE', center);
    % Calculate the derivative of distance squared
    value = (x(1)-x_E(1))*(x(4)-x_E(4)) + (x(2)-x_E(2))*(x(5)-x_E(5)) + (x(3)-x_E(3))*(x(6)-x_E(6));
    isterminal = 1;
    direction = 1;
end

function f = objFunMultipleShoot_parametric_alpha(alpha0, alpha1, TOF, ...
                                                   Delta_v0, alpha_Delta_v0, delta_Delta_v0, msc0, ...
                                                   x2, x3, ...
                                                   betaparameter, Mast, Isp, g0, time_bounds)
    % Objective Function: Maximize distance from Earth at TCA
    
    d2s = 24*60*60;
    
    % Convert parametric times to absolute times
    [t0, t_DSM, t_IMP] = parametric_to_times(alpha0, alpha1, TOF, time_bounds);
    
    % Define the center, frame and bodies
    center = 'Sun';
    frame = 'ECLIPJ2000';
    labels = {'Sun'};
    bodies2 = nbody_init(labels);
    labels_full = {'Sun'; 'Mercury'; 'Venus'; 'Earth'; 'Moon'; ...
                   'Mars Barycenter'; 'Jupiter Barycenter'; 'Saturn Barycenter'; ...
                   'Uranus Barycenter'; 'Neptune Barycenter'; 'Pluto Barycenter'};
    bodies = nbody_init(labels_full);
    
    % Convert parametric launch variables to Cartesian state
    [r1, v1] = parametric_to_cartesian(Delta_v0, alpha_Delta_v0, delta_Delta_v0, t0*d2s, center, frame);
    x1 = [r1; v1];
    
    % Propagate to DSM
    options = odeset('reltol', 1e-12, 'abstol', [1e-6*ones(1,3),1e-9*ones(1,3)]);
    [~, phi1] = ode113(@(t,x) nbody_rhs(t,x,bodies2,frame), [t0 t_DSM]*d2s, x1, options);
    
    % Propagate from DSM to impact
    [~, phi2] = ode113(@(t,x) nbody_rhs(t,x,bodies2,frame), [t_DSM t_IMP]*d2s, x2, options);
    
    % Calculate DSM delta-v
    dv2 = phi2(1,4:6)' - phi1(end,4:6)';
    dv22 = norm(dv2);
    
    % Apply Tsiolkovsky equation
    msc_impact = msc0 * exp(-dv22 / (Isp * g0 * 0.001));
    
    % Calculate dynamic momentum coefficient
    momentumcoeff = msc_impact * betaparameter / Mast;
    
    % Extract asteroid state at impact
    xA_bimp = cspice_spkezr('20099942', t_IMP*d2s, frame, 'NONE', center);
    xA_aimp = [xA_bimp(1:3); xA_bimp(4:6) + x3(4:6)*momentumcoeff];
    
    % Propagate to TCA
    tspan = [t_IMP*d2s, t_IMP*d2s + 7*(365*24*60*60 + 6*60*60)];
    options = odeset('reltol', 1e-12, 'abstol', [1e-6*ones(1,3),1e-9*ones(1,3)],...
        'Events', @(t,x) closest(t,x,center,frame));
    [tt, xx] = ode113(@(t,x) nbody_rhs(t,x,bodies,frame), tspan, xA_aimp, options);
    t_TCA = tt(end);
    xA_TCA = xx(end,:)';
    
    % Calculate distance at TCA
    xE_TCA = cspice_spkezr('Earth', t_TCA, frame, 'NONE', center);
    d = (xA_TCA(1:3) - xE_TCA(1:3));
    f = -(d(1)^2 + d(2)^2 + d(3)^2)^0.5;  % Negative distance (maximize)
end

function [c, ceq] = nonlinconstMultipleShoot_parametric_alpha(alpha0, alpha1, TOF, ...
                                                                Delta_v0, alpha_Delta_v0, delta_Delta_v0, msc0, ...
                                                                x2, x3, ...
                                                                betaparameter, Mast, Isp, g0, time_bounds)
    % Nonlinear constraints with parametric time variables
    
    d2s = 24*60*60;
    
    % Convert parametric times to absolute times
    [t0, t_DSM, t_IMP] = parametric_to_times(alpha0, alpha1, TOF, time_bounds);
    
    % Define the center, frame and bodies
    center = 'Sun';
    frame = 'ECLIPJ2000';
    labels = {'Sun'};
    bodies = nbody_init(labels);
    
    % Convert parametric launch to Cartesian
    [r1, v1] = parametric_to_cartesian(Delta_v0, alpha_Delta_v0, delta_Delta_v0, t0*d2s, center, frame);
    x1 = [r1; v1];
    
    % Initialize equality constraints
    ceq = zeros(12,1);
    
    options = odeset('reltol', 1e-12, 'abstol', [1e-6*ones(1,3),1e-9*ones(1,3)]);
    
    % z1: continuation between flow(x1, t0, t_DSM) and x2 positions
    [~, phi1] = ode113(@(t,x) nbody_rhs(t,x,bodies,frame), [t0 t_DSM]*d2s, x1, options);
    ceq(1:3) = phi1(end,1:3)' - x2(1:3);
    
    % z2: continuation between flow(x2, t_DSM, t_IMP) and x3
    [~, phi2] = ode113(@(t,x) nbody_rhs(t,x,bodies,frame), [t_DSM t_IMP]*d2s, x2, options);
    ceq(4:9) = phi2(end,:)' - x3;
    
    % psi2: match asteroid position at impact
    rvf = cspice_spkezr('20099942', t_IMP*d2s, frame, 'NONE', center);
    ceq(10:12) = x3(1:3) - rvf(1:3);
    
    % Calculate delta-v's
    dv_launch = Delta_v0;
    dv_DSM = norm(phi2(1,4:6)' - phi1(end,4:6)');
    
    % Inequality constraints
    c = [];
    c(end+1) = dv_launch + dv_DSM - 5;          % Total delta-v <= 5 km/s
    c(end+1) = t_IMP - time_bounds.IMPF;        % Impact before window close
    c(end+1) = time_bounds.IMPS - t_IMP;        % Impact after window start
    c(end+1) = t0 - time_bounds.LWC;            % Launch before window close
    c(end+1) = time_bounds.LWO - t0;            % Launch after window open
end

function output = postprocess_parametric_alpha(X_sol, betaparameter, Mast, Isp, g0, time_bounds)
    % Post-processing and visualization
    
    d2s = 24*60*60;
    s2d = 1/d2s;
    center = 'Sun';
    frame = 'ECLIPJ2000';
    labels = {'Sun'; 'Mercury'; 'Venus'; 'Earth'; 'Moon'; ...
              'Mars Barycenter'; 'Jupiter Barycenter'; 'Saturn Barycenter'; ...
              'Uranus Barycenter'; 'Neptune Barycenter'; 'Pluto Barycenter'};
    bodies2 = nbody_init({'Sun'});
    bodies = nbody_init(labels);
    
    % Earth and Moon radii
    Re_v = cspice_bodvrd('EARTH','RADII',3);
    Re = Re_v(1);
    RM_v = cspice_bodvrd('MOON','RADII',3);
    RM = RM_v(1);
    
    % Extract solution variables
    alpha0 = X_sol(1);
    alpha1 = X_sol(2);
    TOF = X_sol(3);
    Delta_v0 = X_sol(4);
    alpha_Delta_v0 = X_sol(5);
    delta_Delta_v0 = X_sol(6);
    msc0 = X_sol(7);
    x2 = X_sol(8:13);
    x3 = X_sol(14:19);
    
    % Convert parametric times to absolute times
    [t0, t_DSM, t_IMP] = parametric_to_times(alpha0, alpha1, TOF, time_bounds);
    
    % Convert parametric launch to Cartesian
    [r1, v1] = parametric_to_cartesian(Delta_v0, alpha_Delta_v0, delta_Delta_v0, t0*d2s, center, frame);
    x1 = [r1; v1];
    
    % Propagate trajectories
    options = odeset('reltol', 1e-12, 'abstol', [1e-6*ones(1,3),1e-9*ones(1,3)]);
    [tt1, phi1] = ode113(@(t,x) nbody_rhs(t,x,bodies2,frame), [t0 t_DSM]*d2s, x1, options);
    [tt2, phi2] = ode113(@(t,x) nbody_rhs(t,x,bodies2,frame), [t_DSM t_IMP]*d2s, x2, options);
    
    % Calculate delta-v's
    rv_earth_launch = cspice_spkezr('Earth', t0*d2s, frame, 'NONE', center);
    dv1 = phi1(1,4:6) - rv_earth_launch(4:6)';
    dv2 = x2(4:6) - phi1(end,4:6)';
    dv11 = norm(dv1);
    dv22 = norm(dv2);
    
    % Mass and momentum calculations
    msc_impact = msc0 * exp(-dv22 / (Isp * g0 * 0.001));
    momentumcoeff = msc_impact * betaparameter / Mast;
    
    % Asteroid trajectory after impact
    xA_bimp = cspice_spkezr('20099942', t_IMP*d2s, frame, 'NONE', center);
    xA_aimp = [xA_bimp(1:3); xA_bimp(4:6) + x3(4:6)*momentumcoeff];
    
    tspan = [t_IMP*d2s, t_IMP*d2s + 10*365*24*60*60];
    options = odeset('reltol', 1e-12, 'abstol', [1e-6*ones(1,3),1e-9*ones(1,3)],...
        'Events', @(t,x) closest(t,x,center,frame));
    [tta, phia] = ode113(@(t,x) nbody_rhs(t,x,bodies,frame), tspan, xA_aimp, options);
    t_TCA = tta(end);
    
    [ttac, phiac] = ode113(@(t,x) nbody_rhs(t,x,bodies,frame), [t_TCA, t_TCA+90*24*3600], phia(end,:)', options);
    phia_m = [phia; phiac];
    
    % Get reference trajectories
    xhA_bimp = cspice_spkezr('20099942', [tt1;tt2]', frame, 'NONE', center)';
    xhE = cspice_spkezr('Earth', [tt1;tt2;tta]', frame, 'NONE', center)';
    
    xA_woImp = cspice_spkezr('20099942', tta', frame, 'NONE', center)';
    xhE2 = cspice_spkezr('Earth', tta', frame, 'NONE', center)';
    xhM = cspice_spkezr('Moon', tta', frame, 'NONE', center)';
    
    % Calculate distances
    dis1 = phia(:,1:3) - xhE2(:,1:3);
    dis2 = xA_woImp(:,1:3) - xhE2(:,1:3);
    d1 = sqrt(sum(dis1.^2, 2));  % With impact
    d2 = sqrt(sum(dis2.^2, 2));  % Without impact
    disM1 = phia(:,1:3) - xhM(:,1:3);
    disM2 = xA_woImp(:,1:3) - xhM(:,1:3);
    dM1 = sqrt(sum(disM1.^2, 2));  % Distance from Moon with Impact
    dM2 = sqrt(sum(disM2.^2, 2));  % Distance from Moon without Impact
  
    
    epoch = datetime('2000-01-01 12:00:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    dateTimes = epoch + seconds(tta');
    
    % FIGURE 1: Distance from Earth at TCA
    figure('Name', 'Distance from Earth at TCA', 'NumberTitle', 'off')
    plot(dateTimes, d1./Re, 'LineWidth', 2, 'DisplayName', 'With Impact')
    hold on
    plot(dateTimes, d2./Re, 'LineWidth', 2, 'DisplayName', 'Without Impact')
    xtickformat('HH:mm:ss');
    xlim([epoch + seconds(tta(end)-3*3600), epoch + seconds(tta(end))])
    legend('Location', 'best')
    ylabel('Distance [R_E]','Interpreter','latex', 'FontSize', 12)
    xlabel('Time [HH:mm:ss]', 'FontSize', 12)
    title('Distance from Earth at TCA - Apophis', 'FontSize', 14)
    grid on
    
    % FIGURE 2: Distance from Moon at TCA
    figure('Name', 'Distance from Moon at TCA', 'NumberTitle', 'off')
    plot(dateTimes, dM1./Re, 'LineWidth', 2, 'DisplayName', 'With Impact')
    hold on
    plot(dateTimes, dM2./Re, 'LineWidth', 2, 'DisplayName', 'Without Impact')
    xtickformat('HH:mm:ss');
    xlim([epoch + seconds(tta(end)-3*3600), epoch + seconds(tta(end))])
    legend('Location', 'best')
    ylabel('Distance [R_E]','Interpreter','latex', 'FontSize', 12)
    xlabel('Time [HH:mm:ss]', 'FontSize', 12)
    title('Distance from Mon at TCA - Apophis', 'FontSize', 14)
    grid on
    
    % FIGURE 3: Mission trajectory in 3D
    figure('Name', '3D Mission Trajectory', 'NumberTitle', 'off')
    plot3(xhE(:,1), xhE(:,2), xhE(:,3), '--', 'LineWidth', 2, 'Color', 'b', 'DisplayName', 'Earth Orbit')
    hold on
    plot3(phi1(:,1), phi1(:,2), phi1(:,3), 'LineWidth', 2, 'DisplayName', 'First Leg (Launch to DSM)')
    plot3(phi2(:,1), phi2(:,2), phi2(:,3), 'LineWidth', 2, 'DisplayName', 'Second Leg (DSM to Impact)')
    plot3(xhA_bimp(length(tt1):end,1), xhA_bimp(length(tt1):end,2), xhA_bimp(length(tt1):end,3), ...
          'LineWidth', 2, 'DisplayName', 'Apophis Before Impact')
    plot3(phia_m(:,1), phia_m(:,2), phia_m(:,3), 'LineWidth', 2, 'Color', 'k', 'DisplayName', 'Apophis After Impact')
    
    % Mark key points
    plot3(phi1(1,1), phi1(1,2), phi1(1,3), 'o', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Departure Point')
    plot3(phi1(end,1), phi1(end,2), phi1(end,3), 'o', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'DSM Point')
    plot3(phi2(end,1), phi2(end,2), phi2(end,3), 'o', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Impact Point')
    plot3(phia(end,1), phia(end,2), phia(end,3), 'o', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Closest Approach')
    
    legend('Location', 'best', 'FontSize', 10)
    grid on
    title('Planetary Defense Mission Overview @Sun-Centered ECLIPJ2000', 'FontSize', 14)
    xlabel('x [km]', 'FontSize', 12)
    ylabel('y [km]', 'FontSize', 12)
    zlabel('z [km]', 'FontSize', 12)
    view(3)
    axis equal
    
    % Print results
    fprintf('\n=== MISSION RESULTS (PARAMETRIC TIME VARIABLES) ===\n');
    output.Launch = cspice_et2utc(t0*d2s,'C',4);
    output.DSM = cspice_et2utc(t_DSM*d2s,'C',4);
    output.IMP = cspice_et2utc(t_IMP*d2s,'C',4);
    output.TCA = cspice_et2utc(t_TCA,'C',4);
    output.alpha0 = alpha0;
    output.alpha1 = alpha1;
    output.TOF = TOF;
    output.Delta_v0 = Delta_v0;
    output.C3 = Delta_v0^2;
    output.alpha_Delta_v0_deg = rad2deg(alpha_Delta_v0);
    output.delta_Delta_v0_deg = rad2deg(delta_Delta_v0);
    output.msc0 = msc0;
    output.msc_impact = msc_impact;
    output.fuel_consumed = msc0 - msc_impact;
    output.momentumcoeff = momentumcoeff;
    output.dvv1 = dv1;
    output.dvv2 = dv2;
    output.dv1_magnitude = dv11;
    output.dv2_magnitude = dv22;
    output.total_dv = dv11 + dv22;
    output.dist = d1(end)/Re;
    
    fprintf('Launch Date: %s\n', output.Launch);
    fprintf('DSM Date: %s\n', output.DSM);
    fprintf('Impact Date: %s\n', output.IMP);
    fprintf('TCA Date: %s\n', output.TCA);
    fprintf('\n--- Parametric Time Variables ---\n');
    fprintf('α₀ (Launch fraction): %.4f\n', alpha0);
    fprintf('α₁ (DSM fraction): %.4f\n', alpha1);
    fprintf('TOF (Total time of flight): %.2f days\n', TOF);
    fprintf('Derived Launch Window Usage: %.1f%% of available window\n', alpha0*100);
    fprintf('Derived DSM Timing: %.1f%% through the mission\n', alpha1*100);
    fprintf('\n--- Launch Parameters ---\n');
    fprintf('Δv₀ (Hyperbolic excess velocity): %.4f km/s\n', Delta_v0);
    fprintf('C3 (Launch Energy): %.4f km²/s²\n', Delta_v0^2);
    fprintf('α_Δv₀ (In-plane angle): %.2f°\n', rad2deg(alpha_Delta_v0));
    fprintf('δ_Δv₀ (Out-of-plane angle): %.2f°\n', rad2deg(delta_Delta_v0));
    fprintf('\n--- Spacecraft Mass Budget ---\n');
    fprintf('Initial Spacecraft Mass (msc0): %.2f kg\n', msc0);
    fprintf('Spacecraft Mass at Impact (msc_impact): %.2f kg\n', msc_impact);
    fprintf('Fuel Consumed: %.2f kg\n', output.fuel_consumed);
    fprintf('Mass Fraction at Impact: %.2f%%\n', (msc_impact/msc0)*100);
    fprintf('Dynamic Momentum Coefficient: %.9f\n', momentumcoeff);
    fprintf('\n--- Delta-V Budget ---\n');
    fprintf('ΔV1 Vector (Launch) [km/s]: [%.6f, %.6f, %.6f]\n', dv1(1), dv1(2), dv1(3));
    fprintf('ΔV1 Magnitude (Launch): %.4f km/s\n', dv11);
    fprintf('ΔV2 Vector (DSM) [km/s]: [%.6f, %.6f, %.6f]\n', dv2(1), dv2(2), dv2(3));
    fprintf('ΔV2 Magnitude (DSM): %.4f km/s\n', dv22);
    fprintf('Total ΔV: %.4f km/s\n', output.total_dv);
    fprintf('ΔV Budget Utilization: %.1f%% of 5.0 km/s limit\n', (output.total_dv/5.0)*100);
    fprintf('\n--- Mission Outcome ---\n');
    fprintf('Distance at TCA: %.4f Earth Radii (%.2f km)\n', output.dist, d1(end));
    fprintf('Improvement over Natural Trajectory: %.2f%%\n', ((d1(end)-d2(end))/d2(end))*100);
    fprintf('===================================================\n');
end


