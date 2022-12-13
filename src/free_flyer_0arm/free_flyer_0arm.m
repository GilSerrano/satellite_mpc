%%  Free Flyer model
%   Simple free flier model taken from Kumar's paper on smooth trajectories
%   for free flyer satellites, e.g. Astrobee.
%   The model paramenters (mass and inertia) were calculated using a CAD
%   software (SolidEdge) for an aluminium cube satellite, with 50cm edge
%   and a 0-joint arm of 50cm.
%   Uses Lie algebra and exponential coordinates to represent orientation.
function model = free_flyer_0arm()

import casadi.*

%% system dimensions
np = 3;
nv = 3;
npee = 3;
nR = 9;
nwb = 3;
nx = np + nv + npee + nR + nwb;
nu1 = 3;
nu2 = 3;
nu = nu1 + nu2;

%% system parameters
M = 79.02;        % mass of satellite [kg]
Jxx= 4.82;  Jyy= 5.40;  Jzz= 5.40;   
Jxy= 0;     Jxz= 0;     Jyz= 0;
J = [Jxx Jxy Jxz; Jxy Jyy Jyz; Jxz Jyz Jzz]; % Inertia of satellite
Jinv = [0.2075, 0, 0;
        0, 0.1852, 0;
        0, 0, 0.1852]; % inv(J)
P_ee = [0.5 0 0]';
    
%% named symbolic variables
p = SX.sym('p', np, 1);         % position of satellite [m]
v = SX.sym('v', nv, 1);         % velocity of satellite [m/s]
pee = SX.sym('pee', np, 1);     % position of end effector [m]
r11 = SX.sym('r11', 1, 1);        % orientation of satellite in rotation matrix
r21 = SX.sym('r21', 1, 1);        % orientation of satellite in rotation matrix
r31 = SX.sym('r31', 1, 1);        % orientation of satellite in rotation matrix
r12 = SX.sym('r12', 1, 1);        % orientation of satellite in rotation matrix
r22 = SX.sym('r22', 1, 1);        % orientation of satellite in rotation matrix
r32 = SX.sym('r32', 1, 1);        % orientation of satellite in rotation matrix
r13 = SX.sym('r13', 1, 1);        % orientation of satellite in rotation matrix
r23 = SX.sym('r23', 1, 1);        % orientation of satellite in rotation matrix
r33 = SX.sym('r33', 1, 1);        % orientation of satellite in rotation matrix
wb = SX.sym('wb', nwb, 1);      % angular velocity of satellite [rad/s]
u1 = SX.sym('u1', nu1, 1);      % force input [N]
u2 = SX.sym('u2', nu2, 1);      % torque input [N/m]

%% (unnamed) symbolic variables
sym_x = vertcat(p, v, pee, r11, r21, r31, r12, r22, r32, r13, r23, r33, wb);
sym_xdot = SX.sym('xdot', nx, 1);
sym_u = vertcat(u1, u2);

%% dynamics
Skew_wb = skew(wb);
r11_dot = r12*wb(3) - r13*wb(2);
r21_dot = r22*wb(3) - r23*wb(2);
r31_dot = r32*wb(3) - r33*wb(2);
r12_dot = -r11*wb(3) + r13*wb(1);
r22_dot = -r21*wb(3) + r23*wb(1);
r32_dot = -r31*wb(3) + r33*wb(1);
r13_dot = r11*wb(2) - r12*wb(1);
r23_dot = r21*wb(2) - r22*wb(1);
r33_dot = r31*wb(2) - r32*wb(1);
R = [r11 r12 r13; r21 r22 r23; r31 r32 r33];

expr_f_expl = vertcat(v, ... % p_dot
                      (1/M) * u1, ... % v_dot
                      v + R * Skew_wb * P_ee, ... % pee_dot
                      r11_dot, ... 
                      r21_dot, ... 
                      r31_dot, ... 
                      r12_dot, ... 
                      r22_dot, ... 
                      r32_dot, ... 
                      r13_dot, ... 
                      r23_dot, ... 
                      r33_dot, ... 
                      Jinv * (u2 - cross(wb, J*wb))); % wb_dot   % simplified explicit dynamics - decoupled


expr_f_impl = expr_f_expl - sym_xdot;   % implicit dynamics

%% constraints
expr_h = sym_u;     % nonlinear constraint h - bounds on input (in this case)

%% cost
W_x = eye(12);
W_u = 0.1*eye(nu);
expr_ext_cost_e = sym_x(7:18)'* W_x * sym_x(7:18); % [pee, r11, r21, r31, r12, r22, r32, r13, r23, r33]
expr_ext_cost = expr_ext_cost_e + sym_u' * W_u * sym_u;
model.W_e = W_x;
model.W = W_x;



%% populate structure
% model.name = 'free_flyer_0arm'; % causing problems in compilation due to file names longer than 63 chars
model.T = 5;     % end time - originally 20
% model.T = 0.9500;
model.nx = nx;
model.nu = nu;
model.sym_x = sym_x;
model.sym_xdot = sym_xdot;
model.sym_u = sym_u;
model.expr_f_expl = expr_f_expl;
model.expr_f_impl = expr_f_impl;
model.expr_h = expr_h;
model.expr_ext_cost = expr_ext_cost;
model.expr_ext_cost_e = expr_ext_cost_e;

end
