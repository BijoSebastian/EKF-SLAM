clc;
clear all;
close all;

%%-------------------------------------------------------------------------
%CONSTANTS
%odometry
m_per_ticks = 0.349;
robot_width = 155.0;
scanner_disp = 30.0;

%Scanner
mounting_angle = -0.06981317007977318;
reduction = 10;

%Kalman Filter
cntrl_mf = 0.35;
cntrl_tf = 0.6;
mea_dist_stddev = 200.0;
mea_ang_stddev = (15.0/180.0)*pi;
d_thresh_match = 500;
sigmax_nl = 100000;
sigmay_nl = 100000;

%Display
xlim = [-500 2500];
ylim = [-500 2500];

figure(1);
hold on;
axis([xlim(1) xlim(2) ylim(1) ylim(2)]);

%%-------------------------------------------------------------------------
%Read file for landmarks.
fileID = fopen('robot_arena_landmarks.txt','r');
formatSpec = '%*s %*s %f %f %*f';
sizeL = [2 inf];
landmarks = fscanf(fileID,formatSpec,sizeL);
fclose(fileID);

%Read file for encoder ticks.
fileID = fopen('robot4_motors.txt','r');
formatSpec = '%*s %*d %d %*d %*d %*d %d %*d %*d %*d %*d %*d %*d %*d';
sizeA = [2 inf];
ticks_A = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
ticks_B = [ticks_A(:,1),ticks_A(:,1:end-1)];
ticks_C = ticks_A-ticks_B;

%Read file for reference trajectory.
fileID = fopen('robot4_reference.txt','r');
formatSpec = '%*s %*d %d %d';
sizeRef = [2 inf];
Ref = fscanf(fileID,formatSpec,sizeRef);
fclose(fileID);
scatter(Ref(1,:),Ref(2,:),10,'filled','green');

%Read file for scan data.
fileID = fopen('robot4_scan.txt','r');
formatSpec = ['%*s %*d %*d' repmat('%d', 1, 660)];
sizeScan = [660 inf];
Scan = fscanf(fileID,formatSpec,sizeScan);
fclose(fileID);

%Plot the map
plot([0,0,2000,2000,0],[0,2000,2000,0,0],'color','black','linewidth',2);
scatter(landmarks(1,:),landmarks(2,:),100,'filled','black');

%%-------------------------------------------------------------------------
%Initial state.
state = [1850.0,1897.0,(213.0/180.0)*pi];
%Initial covariance.
covariance = zeros(3);
%The number of landmarks observed.
num_lndmrks = 0;

%Starting tracking
for i = 1:size(ticks_C,2)

    disp 'time stamp :',i  
    
    %%---------------------------------------------------------------------
    %PREDICTION
    
    %Initialize
    x = state(1);
    y = state(2);
    theta = state(3);
    
    %Ticks to mm
    l = m_per_ticks*ticks_C(1,i);
    r = m_per_ticks*ticks_C(2,i);
    
    %Epsilon Control 
    sigma2_l = (cntrl_mf*l)^2 + (cntrl_tf*(r-l))^2;
    sigma2_r = (cntrl_mf*r)^2 + (cntrl_tf*(r-l))^2;
    eps_cont = diag([sigma2_l, sigma2_r]);
    
    %V matrix
    V = [cos(theta)/2, cos(theta)/2;
         sin(theta)/2, sin(theta)/2;
         -1/robot_width, 1/robot_width];
    %Modify V based on numer of landmarks
    V = [V; zeros(2*num_lndmrks, 2)];
    
    d_c = (l + r)/2.0;
    
    %G matrix
    G = [1, 0, -d_c*sin(theta);
         0, 1, d_c*cos(theta);
         0, 0, 1];
    %Modify G based on numer of landmarks
    G = [G,zeros(3,2*num_lndmrks);
         zeros(2*num_lndmrks,3),eye(2*num_lndmrks)];
    
    %Compute covariance
    covariance = G*covariance*G' + V*eps_cont*V';
    
    %Compute state
    alpha = (r - l)/robot_width;
    
    %weird stuff for scanner disp.
    x = x - scanner_disp*cos(theta);
    y = y - scanner_disp*sin(theta);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    x = x + d_c*cos(theta);
    y = y + d_c*sin(theta);
    theta  = mod((theta + alpha),2*pi);
   
    %weird stuff for scanner disp.
    x = x + scanner_disp*cos(theta);
    y = y + scanner_disp*sin(theta);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    state(1) = x;
    state(2) = y;
    state(3) = theta;
    
    %%---------------------------------------------------------------------
    %%CORRECTION
%%COMMENT START------------------------------------------------------------   
    %The sensing module
    figure(2);
    hold on;
    axis equal;
    
    %Scan data
    len = 1;
    for j = 1:10:size(Scan)
        angle = (j - 330.0) * 0.006135923151543 + mounting_angle; 
        rho(len) = Scan(j,i);
        theta(len) = angle;
        len = len+1;
    end
    [X,Y] = pol2cart(theta',rho');
    h_s(1) = scatter(X,Y,10,'filled','blue');
    
    %Extract lines
    [R, Alpha] = scanner(theta, rho);
    [Xs,Ys] = pol2cart(Alpha, R);
    h_s(2) = scatter(Xs,Ys, 30,'filled','red');
    len_hs = 2;    
    for j = 1:size(R)
        r = R(j);
        alpha = Alpha(j);                                
        vec = [r*cos(alpha),r*sin(alpha)];
        vec = 1000*(vec/norm(vec));
        x1 = -vec(2) + r*cos(alpha);
        y1 = vec(1) + r*sin(alpha);
        x2 = vec(2) + r*cos(alpha);
        y2 = -vec(1) + r*sin(alpha);
        len_hs = len_hs + 1;
        h_s(len_hs) = plot([x2,x1],[y2,y1],'color','green');        
    end
   
    %Correction step
    x = state(1);
    y = state(2);
    theta = state(3);
    for j = 1:size(Xs)  
        if(num_lndmrks ~= 0 && i ~= 1)                            
            Walls = state(4:end);
            r_w = Walls(1:2:end)';
            alpha_w = Walls(2:2:end)';
            
            %Measurement predicition
            r = r_w - (x*cos(alpha_w) +  y*sin(alpha_w));
            alpha  = alpha_w-theta;
            
            %Matching the predicitions with scans. 
            [Xp,Yp] = pol2cart(alpha, r);                                    
            dx = Xs(j) - Xp;
            dy = Ys(j) - Yp;
            dist_match = sqrt((dx.^2) + (dy.^2));
            [small,ind] = min(dist_match);
            if(small < d_thresh_match)
                [alpha_p, r_p] = cart2pol(Xp(ind),Yp(ind));
                [alpha_s, r_s] = cart2pol(Xs(j),Ys(j));
                len_hs = len_hs + 1;
                h_s(len_hs) = plot([Xp(ind),Xs(j)],[Yp(ind),Ys(j)],'color', 'black', 'linewidth', 2);
                
                len_hs = len_hs + 1;
                h_s(len_hs) = scatter(Xp(ind),Yp(ind),30,'filled','black');                

                %Applying corrections.
                H = [-cos(alpha_w(ind)), -sin(alpha_w(ind)), 0.0;
                    0.0, 0.0, -1.0];                    
                %Modify H based on numer of landmarks               
                temp = [1.0, x*sin(alpha_w(ind))+y*cos(alpha_w(ind));
                        0.0, 1.0];
                H = [H, zeros(2,2*(ind - 1)), temp, zeros(2,2*(num_lndmrks - ind))];
                
                %Measurement covariance
                Q = diag([mea_dist_stddev^2, mea_ang_stddev^2]);
                
                %Kalman Gain
                temp = (H*covariance*H') + Q;
                K = (covariance*H')/temp;
                
                %Correcting state
                diff = [r_s - r_p; mod((alpha_s - alpha_p) + pi,2*pi) - pi];%Just to be sure.                                
                state = (state' + K*diff)';
                
                %Correcting covariance                
                covariance = (eye(3+2*num_lndmrks) - K*H)*covariance;
                continue;
            end
        end
        
        %Add the newly observed landmark.
        %Transform to world coordinate system.
        [alpha_s, r_s] = cart2pol(Xs(j),Ys(j));
        alpha_s = alpha_s + theta;
        r_s = r_s + (x*cos(alpha_s) +  y*sin(alpha_s));
        
        %Update state
        state = [state , r_s, alpha_s];
        
        %Update covariance.
        [nr,nc] = size(covariance);
        temp = [sigmax_nl^2, 0.0;
                0.0, sigmay_nl^2];
        covariance = [covariance, zeros(nr,2); zeros(2,nc), temp];
        num_lndmrks = num_lndmrks+1;
        
    end         
%%COMMENT END--------------------------------------------------------------   
    figure(1);
    %Display the covariance ellipse.
    robo_cov = covariance(1:2,1:2);   
    h(1) = drawellipse([state(1), state(2)], robo_cov);
    
    %Display rectified trajectory
    scatter(state(1),state(2),10,'filled','blue');
    
    %Display robot
    x = 150*cos(state(3));
    y = 150*sin(state(3));
    h(2) = plot([state(1),state(1)+x],[state(2),state(2)+y],'color','red','linewidth',2);
    
    alpha = sqrt(covariance(3,3));
    x = [150*cos(mod((state(3)+alpha),2*pi)), 150*cos(mod((state(3)-alpha),2*pi))];
    x = x + state(1);
    y = [150*sin(mod((state(3)+alpha),2*pi)), 150*sin(mod((state(3)-alpha),2*pi))];
    y = y + state(2);
    h(3) = plot([state(1),x,state(1)],[state(2),y,state(2)],'color','red','linewidth',2);
        
    %Display arena.
    len_h = 3;
    for j = 1:num_lndmrks
        r = state(2+(j*2));
        alpha = state(3+(j*2));                        
        vec = [r*cos(alpha),r*sin(alpha)];
        vec = 10000*(vec/norm(vec));
        x1 = -vec(2) + r*cos(alpha);
        y1 = vec(1) + r*sin(alpha);
        x2 = vec(2) + r*cos(alpha);
        y2 = -vec(1) + r*sin(alpha);
        len_h = len_h + 1;
        h(len_h) = plot([x2,x1],[y2,y1],'color','green','linewidth',3);
        %Display the covariance ellipse.
        lnd_mrk_cov = covariance(2+(j*2):3+(j*2),2+(j*2):3+(j*2));
        len_h = len_h + 1;
        h(len_h) = drawellipse(vec, lnd_mrk_cov);        
    end
    
    pause(0.25);
    %Clear display    
%%COMMENT START------------------------------------------------------------
   delete(h_s);
%%COMMENT END--------------------------------------------------------------
    figure(1);
    delete(h);  
    clear h
    clear h_s
end


