%% Drawing the animation of the simulation

pause_t = 0.0001;

satellite_patch_definition;
target_patch_definition;

% Satellite's trajectory
sat_pos = x_sim(1:3,:);
sat_ee_pos = x_sim(7:9,:);
sat_ang = euler;

% Target's trajectory
target_R = R_final;
target_wb = wb_final;
target_pos = (target_R * [0.5 0 0]' + pee_final_2) .*ones(3, N_sim+1);
target_ang = flip(rotm2eul(target_R))' .* ones(3, N_sim+1);
for i = 1:N_sim+1
    target_R = h * target_R * skew(target_wb) + target_R;
    target_ang(:,i) = flip(rotm2eul(target_R));
end

% Path and reference in 3D
figure; hold on;
plot3(sat_ee_pos(1,:), sat_ee_pos(2,:),sat_ee_pos(3,:), '--', 'color', [0.8500 0.3250 0.0980], 'Linewidth', 1.5);
plot3(sat_ee_pos(1,end), sat_ee_pos(2,end),sat_ee_pos(3,end), '*', 'color', [0.8500 0.3250 0.0980], 'Linewidth', 1.5);
xlabel('x'); ylabel('y'); zlabel('z');
grid on;

sat_h = [];
traj_h = [];
target_h = [];
target_traj_h = [];

view_vector = [-1 2 -1];

for k = 1:size(sat_pos,2)
    
    [sat_h, traj_h] = renderSatellite(satV, satF, sat_pos(:,k), sat_ang(:,k), sat_h, traj_h, view_vector);
    [target_h, target_traj_h] = renderTarget(targetV, targetF, target_pos(:,k), target_ang(:,k), target_h, target_traj_h, view_vector);
    
    axis(3.5*[-1, 1, -1, 1, -1, 1]);
    title(sprintf('t = %.3f', ts(k)));
    
    pause(pause_t)
end