% Setup acados environment variables and add folders to path

repo_dir = fileparts(which('setup_acados')); % keep repo in the same parent directory as acados

acados_dir = fullfile(repo_dir, '..', 'acados');
casadi_dir = fullfile(acados_dir, 'external', 'casadi-matlab');
matlab_interface_dir = fullfile(acados_dir, 'interfaces', 'acados_matlab_octave');
mex_template_dir = fullfile(matlab_interface_dir, 'acados_template_mex');

addpath(matlab_interface_dir);
addpath(mex_template_dir);
addpath(casadi_dir);
addpath(genpath('.'))

setenv('ACADOS_INSTALL_DIR', acados_dir);
setenv('ENV_RUN', 'true');