%%  Free Flyer model
%   Simple free flier model taken from Kumar's paper on smooth trajectories
%   for free flyer satellites, e.g. Astrobee.
%   The model paramenters (mass and inertia) were calculated using a CAD
%   software (SolidEdge) for an aluminium cube satellite, with 50cm edge
%   and a 0-joint arm of 50cm.
%   Uses Lie algebra and exponential coordinates to represent orientation.
function model = free_flyer_simple()

import casadi.*

%% system dimensions
np = 3;
nv = 3;
nksi = 3;
nwb = 3;
nx = np + nv + nksi + nwb;
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
    
%% named symbolic variables
p = SX.sym('p', np, 1);         % position of satellite [m]
ksi = SX.sym('ksi', nksi, 1);   % orientation of satellite in exponential coordinates
v = SX.sym('v', nv, 1);         % velocity of satellite [m/s]
wb = SX.sym('wb', nwb, 1);      % angular velocity of satellite [rad/s]
u1 = SX.sym('u1', nu1, 1);              % force input [N]
u2 = SX.sym('u2', nu2, 1);              % torque input [N/m]

%% (unnamed) symbolic variables
sym_x = vertcat(p, v, ksi, wb);
sym_xdot = SX.sym('xdot', nx, 1);
sym_u = vertcat(u1, u2);

%% dynamics
ksi_norm = norm(ksi) + eps;     % Norm of the exponential coordinates (added small epsilon to prevent sigularity ksi = [0 0 0])
Skew_ksi = skew(ksi);           % Skew-symmetric matrix of the exponential coordinates

Jr_inv = eye(3) + 0.5*Skew_ksi + (1/ksi_norm^2 - (1+cos(ksi_norm))/(2*ksi_norm*sin(ksi_norm)) )*Skew_ksi^2;     % Right Jacobian inverse

% expr_f_expl = vertcat(v, ...
%                       u1, ...
%                       Jr_inv * wb, ...
%                       u2); % simplified explicit dynamics - double integrator
expr_f_expl = vertcat(v, ...
                      (1/M) * u1, ...
                      Jr_inv * wb, ...
                      Jinv * (u2 - cross(wb, J*wb)));   % simplified explicit dynamics - decoupled


expr_f_impl = expr_f_expl - sym_xdot;   % implicit dynamics

%% constraints
expr_h = sym_u;     % nonlinear constraint h - bounds on input (in this case)

%% cost
W_x = eye(6); %diag([1e3, 1e3, 1e-2, 1e-2]);
% W_x(7:12,7:12) = 0.1*eye(6);
W_u = 0.1*eye(nu); %1e-2;
expr_ext_cost_e = [sym_x(1:3); sym_x(7:9)]'* W_x * [sym_x(1:3); sym_x(7:9)];
expr_ext_cost = expr_ext_cost_e + sym_u' * W_u * sym_u;
model.cost_expr_ext_cost = expr_ext_cost;
model.cost_expr_ext_cost_e = expr_ext_cost_e;
% nonlinear least sqares
cost_expr_y = [eye(3) zeros(3,9); zeros(3,6) eye(3) zeros(3)]*sym_x;
model.cost_expr_y_e = [eye(3) zeros(3,9); zeros(3,6) eye(3) zeros(3)]*sym_x;
model.W_e = W_x;
W = W_x;
% cost_expr_y = vertcat(sym_x, sym_u);
% W = blkdiag(W_x, W_u);
% model.cost_expr_y_e = sym_x;
% model.W_e = W_x;



%% populate structure
model.name = 'free_flyer';
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

model.cost_expr_y = cost_expr_y;
model.W = W;

end
