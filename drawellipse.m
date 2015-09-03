function [h] = drawellipse(mu,S)
    %# This functions draw an ellipse
    [U,D] = eig(S); % U = eigenvectors, D= diagonal matrix of eigenvalues.
    if(D(1,1)>D(2,2))
        main_axis_angle = atan2(U(2,1),U(1,1));
        radius1 = sqrt(D(1,1));
        radius2 = sqrt(D(2,2));        
    else
        main_axis_angle = atan2(U(2,2),U(1,2));
        radius1 = sqrt(D(2,2));
        radius2 = sqrt(D(1,1)); 
    end
    
    start_angle = 0.0;
    end_angle = 2 * pi;
    x =  [];
    y = [];
    ax = radius1 * cos(main_axis_angle);
    ay = radius1 * sin(main_axis_angle);
    bx = - radius2 * sin(main_axis_angle);
    by = radius2 * cos(main_axis_angle);
    N_full = 40;
%     N = int(ceil((end_angle - start_angle) / (2 * pi) * N_full));
%     N = max(N, 1);
%     increment = (end_angle - start_angle) / N;
increment = (end_angle - start_angle) / N_full;
    
    for i =1:(N_full + 1)
        a = start_angle + i * increment;
        c = cos(a);
        s = sin(a);
        x = [x,c*ax + s*bx + mu(1)];
        y = [y,- c*ay - s*by + mu(2)];
    end
    h = plot(x,y);
      

end
