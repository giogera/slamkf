clear all
close all
clc

load 'lidar_scan.mat'

%% Sampling times and estimates of covariance matrices Q, Q_turn, R

% Vector of sampling times
Ts = diff(TimeStamp);

% Covariance matrices estimation
Q = cov([Uf(373:674)' Ua(373:674)']); % second hallway
Qturn = cov([Uf(1261:1427)' Ua(1261:1427)']); % fourth curve
R = [8e-2, 0.0001; 0.0001, 3e-2];

% Make Q and Qturn white
Q(1,2) = 0; Q(2,1) = 0; 
Qturn(1,2) = 0; Qturn(2,1) = 0;

wturn = 0.1897;

%% Initialization
P_odom = (1/1000 * P_odom)';   % conversion into meters
N = length(P_odom);
N_meas = size(noisyRangeData,1); % number of measurements for each scan

% Preallocate the vectors for estimates storing
PoseEst = zeros(N,3); % pose of the robot over time
LandmarksEst =[]; % position of the landmarks

% State and Covariance Initialization
Z0 = P_odom(1,:)';
P0 = diag([0.01 0.01 0.0001]);

% Thresholds
tau1 = 5.9915;
tau2 = 15;

eta = 10;

%% Extended Kalman Filter

t = 1;
lid = 0; % number of identified landmarks
nla(t) = lid;

% Preallocations
Hl = zeros(2,3); % landmark measurement Jacobian
hl = zeros(2,1); % landmark expected position
Kl = zeros(2*lid+3,2); % landmark Kalman gain
Zc = zeros(3,1); % correction step state
Zp = zeros(3,1); % prediction step state

P = P0;
Zc = Z0;

% Map coordinates initialization
Xmappose = zeros(N_meas,N);
Ymappose = zeros(N_meas,N);

% EKF Iteration

figure
hold on

for t = 1:N 

    % --- Correction Step --- %

    range = flip(noisyRangeData(:,t));
    [m1, m2, Xc, Yc] = LandmarkSeen(t,range,Zc,ang_span,N_meas);

    % if t==400
    %     figure
    %     hold on
    %     plot(Xmappose,Ymappose,'.r')
    %     %plot(Xmapodom,Ymapodom,'.g')
    %     hold off
    %     keyboard
    % end

    if lid == 0 % no landmarks have been seen

        for j = 1:size(m1,1)

            m = [m1(j);m2(j)];
            lid = lid + 1;
            nla(t) = lid;

            Zc = [Zc; Zc(1) + m(1)*cos(m(2) + Zc(3)); Zc(2) + m(1)*sin(m(2) + Zc(3))];
            P = blkdiag(P, eta*eye(2)); %augmented covariance matrix with the new landmark

            dx = Zc(2*lid+2) - Zc(1);
            dy = Zc(2*lid+3) - Zc(2);
            den = dx^2 + dy^2;

            % Landmark Expected Position (range and angle)
            hl(1) = sqrt(den);
            hl(2) = atan2(dy, dx) - Zc(3);

            % Measurement Jacobian
            Hl = zeros(2,2*lid+3);
            Hl(1,1:2) = [-dx/sqrt(den) -dy/sqrt(den)];
            Hl(1,(2*lid+2):(2*lid+3)) = [dx/sqrt(den) dy/sqrt(den)];
            Hl(2,1:3) = [dy/den -dx/den -1];
            Hl(2,(2*lid+2):(2*lid+3)) = [-dy/den dx/den];

            % Kalman Gain
            Kl = P*Hl'*inv(Hl*P*Hl' + R);

            % Correct the predicted State Estimate and the Covariance Matrix
            Zc = Zc + Kl*[m(1)-hl(1); wrapToPi(m(2)-hl(2))];
            P = P - P*Hl'*Kl';

        end

    else

        % do the correction when I have at least one landmark in memory

        lt = length(m1);
        for j = 1:lt

            m = [m1(j) m2(j)]'; % we assume this measurement is associated to landmark Lk
            distance_vector = zeros(1,lid);

            for k = 1:lid

                arg1 = (Zc(2*k+2)-Zc(1));
                arg2 = (Zc(2*k+3)-Zc(2));
                deltajk = [m(1)-sqrt(arg1^2+arg2^2);
                    wrapToPi(m(2)-atan2(arg2,arg1)+Zc(3))];

                dx = Zc(2*k+2) - Zc(1);
                dy = Zc(2*k+3) - Zc(2);
                den = dx^2 + dy^2;

                Hk = zeros(2,2*lid+3);
                Hk(1,1:2) = [-dx/sqrt(den) -dy/sqrt(den)];
                Hk(1,(2*k+2):(2*k+3)) = [dx/sqrt(den) dy/sqrt(den)];
                Hk(2,1:3) = [dy/den -dx/den -1];
                Hk(2,(2*k+2):(2*k+3)) = [-dy/den dx/den];

                distance = (deltajk)'* inv((Hk)*P*(Hk)'+R) * (deltajk);

                distance_vector(k) = distance;

            end

            [dj_star_val,dj_star_pos] = min(distance_vector); % it gives back the position and value of the minimum distance landmark

            if dj_star_val < tau1

                dx = Zc(2*dj_star_pos+2) - Zc(1);
                dy = Zc(2*dj_star_pos+3) - Zc(2);
                den = dx^2 + dy^2;

                % Landmark Expected Position (range and angle)
                hl(1) = sqrt(den);
                hl(2) = atan2(dy, dx) - Zc(3);

                % Measurement Jacobian
                Hl = zeros(2,2*lid+3);
                Hl(1,1:2) = [-dx/sqrt(den) -dy/sqrt(den)];
                Hl(1,(2*dj_star_pos+2):(2*dj_star_pos+3)) = [dx/sqrt(den) dy/sqrt(den)];
                Hl(2,1:3) = [dy/den -dx/den -1];
                Hl(2,(2*dj_star_pos+2):(2*dj_star_pos+3)) = [-dy/den dx/den];

                % Kalman Gain
                Kl = P*Hl'*inv(Hl*P*Hl' + R);

                % Correct the predicted State Estimate and the Covariance Matrix
                Zc = Zc + Kl*[m(1)-hl(1); wrapToPi(m(2)-hl(2))];
                P = P - P*Hl'*Kl';

            end

            if dj_star_val > tau2  % new landmark found

                lid = lid + 1;
                Zc = [Zc; Zc(1) + m(1)*cos(m(2) + Zc(3)); Zc(2) + m(1)*sin(m(2) + Zc(3))];
                P = blkdiag(P, eta*eye(2)); % augmented covariance matrix with the new landmark

                dx = Zc(2*lid+2) - Zc(1);
                dy = Zc(2*lid+3) - Zc(2);
                den = dx^2 + dy^2;

                % Landmark Expected Position (range and angle)
                hl(1) = sqrt(den);
                hl(2) = atan2(dy, dx) - Zc(3);

                % Measurement Jacobian
                Hl = zeros(2,2*lid+3);
                Hl(1,1:2) = [-dx/sqrt(den) -dy/sqrt(den)];
                Hl(1,(2*lid+2):(2*lid+3)) = [dx/sqrt(den) dy/sqrt(den)];
                Hl(2,1:3) = [dy/den -dx/den -1];
                Hl(2,(2*lid+2):(2*lid+3)) = [-dy/den dx/den];

                % Kalman Gain
                Kl = P*Hl'*inv(Hl*P*Hl' + R);

                % Correct the predicted State Estimate and the Covariance Matrix
                Zc = Zc + Kl*[m(1)-hl(1); wrapToPi(m(2)-hl(2))];
                P = P - P*Hl'*Kl';

            end
        end

    end

    % Store the corrected State Estimate and Landmark estimate
    PoseEst(t,:) = Zc(1:3)';
    for j=1:lid
        LandmarksEst(j,1) = Zc(2*j+2);
        LandmarksEst(j,2) = Zc(2*j+3);
    end

    % --- Prediction Step --- %

    if t == N
        break;
    end

    Zp = Zc;

    % Position update
    Zp(1) = Zc(1) + Ts(t)*Uf(t)*cos(Zc(3));
    Zp(2) = Zc(2) + Ts(t)*Uf(t)*sin(Zc(3));
    Zp(3) = Zc(3) + Ts(t)*Ua(t);

    % Covariance of the process disturbance
    if abs(Ua(t)) > wturn
        Qp = Qturn;
    else
        Qp = Q;
    end

    % Robot Pose Jacobians
    Fu = [1 0 -Ts(t)*Uf(t)*sin(Zc(3)); 0 1 Ts(t)*Uf(t)*cos(Zc(3)); 0 0 1];
    Gu = [Ts(t)*cos(Zc(3)) 0; Ts(t)*sin(Zc(3)) 0; 0 Ts(t)];

    % "Complete" Jacobians
    F = blkdiag(Fu, eye(2*lid)); 
    G = [Gu; zeros(2*lid,2)]; 
    clear Fu Gu

    P = F*P*F' + G*Qp*G';
    Zc = Zp;

    plot(PoseEst(t,1),PoseEst(t,2),'go')
    % plot(P_odom(t,1),P_odom(t,2),'mo')

    if t==1
        for nl=1:lid
            handl(nl)=plot(LandmarksEst(nl,1),LandmarksEst(nl,2),'*c');
        end
    else
        for nl=1:nla(t-1)
            set(handl(nl),'XData',LandmarksEst(nl,1),'YData',LandmarksEst(nl,2));
        end
        for nl=nla(t-1)+1:lid
            handl(nl)=plot(LandmarksEst(:,1),LandmarksEst(:,2),'*c');
        end
    end

    pause(0.001)
    nla(t)=lid;

    [x,y] = scanmap(range,Zc,ang_span,N_meas);
    Xmappose(1:length(x),t) = x';
    Ymappose(1:length(x),t) = y';

    if mod(t-1,10) == 0
        plot(x,y,'.r')
    end
    
end

figure
plot(Xmappose,Ymappose, '.r')
hold on
for i=1:length(LandmarksEst)
    plot(LandmarksEst(i,1),LandmarksEst(i,2),'.c')
end
 
%% Functions

function [yy, alpha, Xc, Yc] = LandmarkSeen(t,x,Zc,ang_span,N_meas)

% Parameters for landmark selection
ValProminence = 0.035;
threshold_left = -0.01;
threshold_right = 0.01;

ra = ang_span/(N_meas - 1); 
spaz = -ang_span/2:ra:ang_span/2; % sign changed to flip the angles on the kalman filter 
tx = 1:1:N_meas;
[TF,P] = islocalmin(x,'MinProminence',ValProminence);

% Filter to avoid values equal to 0 (measurement too distant)
z=x(TF);
temp=tx(TF);
y=[];
ty=[];
for i=1:length(z)
    if z(i)>0        
        y = [y; z(i)];
        ty = [ty; temp(i)];
    end    
end

% Filter on the slope

yy = [];         % distance associated to the minimum (meas.range)
tyy = [];        % i-th scan associated to the local minimum yy considered
alpha = [];      % angular corresponding of tyy
threshold_list_left = [];
threshold_list_right = [];

for i = 1:length(y)

    if ty(i) > 3 && ty(i)< N_meas-2
        m_left1 = (y(i)-x(ty(i)-1))/(ty(i)-tx(ty(i)-1));
        m_left2 = (y(i)-x(ty(i)-2))/(ty(i)-tx(ty(i)-2));
        m_left = (m_left1+m_left2)/2;

        m_right1 = (x(ty(i)+1)-y(i))/(tx(ty(i)+1)-ty(i));
        m_right2 = (x(ty(i)+2)-y(i))/(tx(ty(i)+2)-ty(i));
        m_right = (m_right1+m_right2)/2;

        threshold_list_left = [threshold_list_left;m_left];
        threshold_list_right = [threshold_list_right;m_right];

        if m_left < threshold_left && m_right > threshold_right
            yy = [yy; y(i)];
            tyy = [tyy; ty(i)];
            alpha = [alpha; spaz(ty(i))];
        end

    end
end

% Compute the coordinates of the landmarks
Xc = zeros(length(yy),1);
Yc = zeros(length(yy),1);
for j=1:length(yy)
    Xc(j) = Zc(1)+yy(j)*cos(spaz(tyy(j))+Zc(3));
    Yc(j) = Zc(2)+yy(j)*sin(spaz(tyy(j))+Zc(3));
end

end

function [Xc, Yc] = scanmap(x,Z,ang_span,N_meas)

ra = ang_span/(N_meas - 1);
spaz = -ang_span/2:ra:ang_span/2; % sign changed to flip the angles on the kalman filter 
tx = 1:1:N_meas;

% Filter to avoid values equal to 0 (measurement too distant)
y=[];
ty=[];
for i=1:length(x)
    if x(i)>0        
        y = [y; x(i)];
        ty = [ty; tx(i)];
    end    
end

Xc = zeros(length(y),1);
Yc = zeros(length(y),1);
for j=1:length(y)
    Xc(j) = Z(1)+y(j)*cos(spaz(ty(j))+Z(3));
    Yc(j) = Z(2)+y(j)*sin(spaz(ty(j))+Z(3));
end

end