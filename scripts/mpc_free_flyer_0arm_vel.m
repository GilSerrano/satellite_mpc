%% Closed loop MPC of simple free flyer
clear; clc; 

check_acados_requirements()

% initial state
x0 = [0; 0; 0; 0; 0; 0; 0.5; 0; 0; 1; 0; 0; 0; 1; 0; 0; 0; 1; 0; 0; 0];  % start at stable position

% reference
pee_final = [0 1 0]';
wb_x_ref = 1; % rotation in the end can only be around the x-axis (axis with arm)
wb_final = wb_x_ref*[1 0 0]';
% 
% yref = wb_final;
% yref_e = wb_final;

yref = zeros(6, 1);
yref_e = zeros(6, 1);

yref(1:3) = pee_final;
yref(4:6) = wb_final;
yref_e(1:3) = pee_final;
yref_e(4:6) = wb_final;

%% discretization
h = 0.1; % sampling time = length of first shooting interval
N = 20; % number of shooting intervals
% nonuniform discretization
shooting_nodes = [0.0 h, 5*h*(1:N-1)];
% shooting_nodes = [0.0 0.01, 0.05*(1:N-1)];
T = shooting_nodes(end);

nlp_solver = 'sqp'; % sqp, sqp_rti
qp_solver = 'partial_condensing_hpipm';
% full_condensing_hpipm, partial_condensing_hpipm, full_condensing_qpoases, full_condensing_daqp
qp_solver_cond_N = 5; % for partial condensing

% we add some model-plant mismatch by choosing different integration
% methods for model (within the OCP) and plant:

% integrator model
model_sim_method = 'erk';
model_sim_method_num_stages = 1;
model_sim_method_num_steps = 2;

% integrator plant
plant_sim_method = 'irk';
plant_sim_method_num_stages = 3;
plant_sim_method_num_steps = 3;

%% model dynamics
model = free_flyer_0arm;
nx = model.nx;
nu = model.nu;

%% model to create the solver
ocp_model = acados_ocp_model();

%% acados ocp model
ocp_model.set('T', T);
% symbolics
ocp_model.set('sym_x', model.sym_x);
ocp_model.set('sym_u', model.sym_u);
ocp_model.set('sym_xdot', model.sym_xdot);

% external cost function
ocp_model.set('cost_type', 'ext_cost');
ocp_model.set('cost_type_e', 'ext_cost');

costQ = [1*eye(3), zeros(3); zeros(3) 10*eye(3)];
costQv = 100*eye(3);
costR = 1*eye(6);
expr_ext_cost_e = ([ocp_model.model_struct.sym_x(7:9); ocp_model.model_struct.sym_x(19:21)] - yref)'* costQ * ([ocp_model.model_struct.sym_x(7:9); ocp_model.model_struct.sym_x(19:21)] - yref) + (ocp_model.model_struct.sym_x(4:6))' * costQv * (ocp_model.model_struct.sym_x(4:6));
expr_ext_cost = expr_ext_cost_e + ocp_model.model_struct.sym_u' * costR * ocp_model.model_struct.sym_u;
ocp_model.set('cost_expr_ext_cost', expr_ext_cost);
ocp_model.set('cost_expr_ext_cost_e', expr_ext_cost_e);

% dynamics
ocp_model.set('dyn_type', 'explicit');
ocp_model.set('dyn_expr_f', model.expr_f_expl);

% constraints
ocp_model.set('constr_type', 'auto');
ocp_model.set('constr_expr_h', model.expr_h);
U_max = 10*ones(nu,1);
ocp_model.set('constr_lh', -U_max); % lower bound on h
ocp_model.set('constr_uh', U_max);  % upper bound on h
ocp_model.set('constr_x0', x0);

%% acados ocp set opts
ocp_opts = acados_ocp_opts();
ocp_opts.set('param_scheme_N', N);
ocp_opts.set('shooting_nodes', shooting_nodes);

ocp_opts.set('nlp_solver', nlp_solver);
ocp_opts.set('sim_method', model_sim_method);
ocp_opts.set('sim_method_num_stages', model_sim_method_num_stages);
ocp_opts.set('sim_method_num_steps', model_sim_method_num_steps);

ocp_opts.set('qp_solver', qp_solver);
ocp_opts.set('qp_solver_cond_N', qp_solver_cond_N);
% ... see ocp_opts.opts_struct to see what other fields can be set

ocp_opts.set('print_level', 1);

%% create ocp solver
ocp = acados_ocp(ocp_model, ocp_opts);

x_traj_init = zeros(nx, N+1);
u_traj_init = zeros(nu, N);


%% plant: create acados integrator
% acados sim model
sim_model = acados_sim_model();
sim_model.set('name', 'plant');
sim_model.set('T', h);

sim_model.set('sym_x', model.sym_x);
sim_model.set('sym_u', model.sym_u);
sim_model.set('sym_xdot', model.sym_xdot);
sim_model.set('dyn_type', 'implicit');
sim_model.set('dyn_expr_f', model.expr_f_impl);

% acados sim opts
sim_opts = acados_sim_opts();
sim_opts.set('method', plant_sim_method);
sim_opts.set('num_stages', plant_sim_method_num_stages);
sim_opts.set('num_steps', plant_sim_method_num_steps);

sim = acados_sim(sim_model, sim_opts);

%% Simulation
t_sim = 100;
N_sim = t_sim/h;
x0 = [0; 0; 0; 0; 0; 0; 0.5; 0; 0; 1; 0; 0; 0; 1; 0; 0; 0; 1; 0; 0; 0];  % start at stable position

x_sim = zeros(nx, N_sim+1);
u_sim = zeros(nu, N_sim);

x_sim(:,1) = x0;

for i=1:N_sim
    % update initial state
    x0 = x_sim(:,i);
    ocp.set('constr_x0', x0);

    % solve
    ocp.solve();

    % get solution
    u0 = ocp.get('u', 0);
    status = ocp.get('status'); % 0 - success
    if status 
        disp('not successful');
    end

    % set initial state
    sim.set('x', x0);
    sim.set('u', u0);

    % solve
    sim_status = sim.solve();
    if sim_status ~= 0
        disp(['acados integrator returned error status ', num2str(sim_status)])
    end

    % get simulated state
    x_sim(:,i+1) = sim.get('xn');
    u_sim(:,i) = u0;
end


%% Plots
ts = linspace(0, N_sim*h, N_sim+1);
States = {'p', 'v', 'pee', 'r1', 'r2', 'r3', 'wb'};

% com position
figure;
plot(ts, x_sim(1,:), 'color', [0 0.4470 0.7410]);hold on;
plot(ts, x_sim(2,:), 'color', [0.8500 0.3250 0.0980]);
plot(ts, x_sim(3,:), 'color', [0.9290 0.6940 0.1250]);
ylabel('CoM position'); xlabel('t [s]');
legend('p_x', 'p_y', 'p_z', 'location', 'best');
grid on; 

% com velocity
figure;
plot(ts, x_sim(4,:), 'color', [0 0.4470 0.7410]);hold on;
plot(ts, x_sim(5,:), 'color', [0.8500 0.3250 0.0980]);
plot(ts, x_sim(6,:), 'color', [0.9290 0.6940 0.1250]);
ylabel('CoM velocity'); xlabel('t [s]');
legend('v_x', 'v_y', 'v_z', 'location', 'best');
grid on; 

% EE position
figure;
plot(ts, x_sim(7,:), 'color', [0 0.4470 0.7410]);hold on;
plot(ts, x_sim(8,:), 'color', [0.8500 0.3250 0.0980]);
plot(ts, x_sim(9,:), 'color', [0.9290 0.6940 0.1250]);
plot(ts, yref(1)*ones(length(ts),1), '--', 'color', [0 0.4470 0.7410]);
plot(ts, yref(2)*ones(length(ts),1), '--', 'color', [0.8500 0.3250 0.0980]);
plot(ts, yref(3)*ones(length(ts),1), '--', 'color', [0.9290 0.6940 0.1250]);
ylabel('EE position'); xlabel('t [s]');
legend('pee_x', 'pee_y', 'pee_z', 'location', 'best');
grid on; 

% orientation
R_sim = zeros(3,3,N_sim+1);
euler = zeros(3,N_sim+1);
determ = zeros(N_sim+1, 1);
for i = 1:N_sim+1
    R_sim(1,1,i) = x_sim(10,i);
    R_sim(2,1,i) = x_sim(11,i);
    R_sim(3,1,i) = x_sim(12,i);
    R_sim(1,2,i) = x_sim(13,i);
    R_sim(2,2,i) = x_sim(14,i);
    R_sim(3,2,i) = x_sim(15,i);
    R_sim(1,3,i) = x_sim(16,i);
    R_sim(2,3,i) = x_sim(17,i);
    R_sim(3,3,i) = x_sim(18,i);
    determ(i) = det(R_sim(:,:,i));
    euler_aux = rotm2eul(R_sim(:,:,i)); % [yaw pitch roll]
    euler(:,i) = flip(euler_aux); % [roll pitch yaw]
end

% Rotation reference
figure;
plot(ts, euler(1,:), 'color', [0 0.4470 0.7410]);hold on;
plot(ts, euler(2,:), 'color', [0.8500 0.3250 0.0980]);
plot(ts, euler(3,:), 'color', [0.9290 0.6940 0.1250]);
ylabel('Orientation'); xlabel('t [s]');
legend('roll', 'pitch', 'yaw', 'location', 'best');
grid on;

% Angular velocity
figure;
plot(ts, x_sim(19,:), 'color', [0 0.4470 0.7410]);hold on;
plot(ts, x_sim(20,:), 'color', [0.8500 0.3250 0.0980]);
plot(ts, x_sim(21,:), 'color', [0.9290 0.6940 0.1250]);
plot(ts, yref(4)*ones(length(ts),1), '--', 'color', [0 0.4470 0.7410]);
plot(ts, yref(5)*ones(length(ts),1), '--', 'color', [0.8500 0.3250 0.0980]);
plot(ts, yref(6)*ones(length(ts),1), '--', 'color', [0.9290 0.6940 0.1250]);
ylabel('Angular velocity'); xlabel('t [s]');
legend('wb_x', 'wb_y', 'wb_z', 'location', 'best');
grid on;

u_vec = [u_sim'; u_sim(:,end)'];
figure; hold on;
stairs(ts, u_vec(:,1));
stairs(ts, u_vec(:,2));
stairs(ts, u_vec(:,3));
ylabel('u1'); xlabel('t [s]');
legend('u1_x', 'u1_y', 'u1_z', 'location', 'best');
grid on;
figure; hold on;
stairs(ts, u_vec(:,4));
stairs(ts, u_vec(:,5));
stairs(ts, u_vec(:,6));
ylabel('u2'); xlabel('t [s]');
legend('u2_x', 'u2_y', 'u2_z', 'location', 'best');
grid on;
