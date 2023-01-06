function [sat_h,traj_h] = renderSatellite(V,F,p,ang,sat_h,traj_h, view_vector)
                                 
    V = rotate(V,ang);
    V = translate(V,p');

    if isempty(sat_h) && isempty(traj_h)
    figure(1)
    traj_h = plot3(p(1),p(2),p(3),'--','color', [0 0.4470 0.7410],'Linewidth',1.5); 
    hold on;
    sat_h = patch('Vertices', V', 'Faces', F, 'FaceColor','#808080');
%     view(3);
    view(view_vector);
    
    else
        set(sat_h,'Vertices',V','Faces',F);
        set(traj_h,'Xdata',[get(traj_h,'Xdata'),p(1)]);
        set(traj_h,'Ydata',[get(traj_h,'Ydata'),p(2)]);
        set(traj_h,'Zdata',[get(traj_h,'Zdata'),p(3)]);
    end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function V = rotate(V,ang)

    phi = ang(1); %roll
    theta = ang(2); % pitch
    psi = ang(3); % yaw
    
    R_roll = [...
          1, 0, 0;...
          0, cos(phi), -sin(phi);...
          0, sin(phi), cos(phi)];
    R_pitch = [...
          cos(theta), 0, sin(theta);...
          0, 1, 0;...
          -sin(theta), 0, cos(theta)];
    R_yaw = [...
          cos(psi), -sin(psi), 0;...
          sin(psi), cos(psi), 0;...
          0, 0, 1];
    
    R = R_yaw*R_pitch*R_roll;
    V = R*V;
%     V = V';
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function V = translate(V,p)

  V = V + repmat(p',1,size(V,2));
  
end


