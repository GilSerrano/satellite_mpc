%% Drawing the animation of the simulation

pause_t = 0.0001;

satellite_patch_definition;

sat_pos = x_sim(1:3,:);
sat_ee_pos = x_sim(7:9,:);
sat_ang = euler;

% Path and reference in 3D
figure; hold on;
plot3(sat_ee_pos(1,:), sat_ee_pos(2,:),sat_ee_pos(3,:), '--', 'color', [0.8500 0.3250 0.0980], 'Linewidth', 1.5);
plot3(sat_ee_pos(1,end), sat_ee_pos(2,end),sat_ee_pos(3,end), '*', 'color', [0.8500 0.3250 0.0980], 'Linewidth', 1.5);
xlabel('x'); ylabel('y'); zlabel('z');
grid on;

sat_h = [];
traj_h = [];
for k = 1:size(sat_pos,2)
    
    [sat_h, traj_h] = renderSatellite(satV, satF, sat_pos(:,k), sat_ang(:,k),sat_h,traj_h);
    
    axis(2*[-1.5, 1.5, -1.5, 1.5, -1.5, 1.5]);
%     axis([sat_pos(1,k)-1.5, sat_pos(1,k)+1.5, sat_pos(2,k)-0.5, sat_pos(2,k)+0.5, sat_pos(3,k)-1.5, sat_pos(3,k)+1.5]);
    title(sprintf('t = %.3f', ts(k)));
    
    pause(pause_t)
end