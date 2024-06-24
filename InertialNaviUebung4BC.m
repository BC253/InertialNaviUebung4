clear all
% close all
clc
addpath("./Funktion")

w_E = [0;0;7.292115e-5];
Omega_iee = [0 -w_E(3) 0
             w_E(3) 0 0
             0 0 0];
%% startwert

% 读取IMU数据 ax ay az(m/s^2) gx gy gz(rad/s) 0.02s
IMU.raw = csvread('IMU.csv', 1, 0);
t_imu = IMU.raw(:,1);
IMUdata(:,1:3) = IMU.raw(:,2:4);
IMUdata(:,4:6) = IMU.raw(:,5:7);
%GNSS 数据
GNSS.raw = csvread('PVA-GNSS.csv', 1, 0);% 1s
gnssp_e = GNSS.raw(:,2:4);  %%初始数据 startwert
gnssv_e = GNSS.raw(:,8:10);%%初始数据 startwert
gnssv_n = GNSS.raw(:,11:13);
gnsslla = GNSS.raw(:,5:7);
gnsslla = deg2rad(gnsslla(:,1:2));
gnssRPY = deg2rad(GNSS.raw(:,14:16));
t_gnss = GNSS.raw(:,1);
% REF数据
REF.raw = csvread('PVA-REF.csv', 1, 0);%0.2s
refp = REF.raw(:,2:4);
refv = REF.raw(:,8:10);
refv_n = REF.raw(:,11:13);
reft = REF.raw(:,1);

%DCM t0 berechnen
Cne_0 = C(3,-gnsslla(1,2))*C(2,gnsslla(1,1)+pi/2);
Cbn_0 = C(3,-gnssRPY(1,3))*C(2,-gnssRPY(1,2))*C(1,-gnssRPY(1,1));
Cpe_0 = Cne_0*Cbn_0;
Omega_iep = inv(Cpe_0)*Omega_iee*Cpe_0;
w_ieP(1,1) = Omega_iep(3,2);
w_ieP(2,1) = Omega_iep(1,3);
w_ieP(3,1) = Omega_iep(2,1);
                                          

%DCM t1 berechnen
Cne_1 = C(3,-gnssp_e(2,2))*C(2,gnssp_e(2,1)+pi/2);
Cbn_1 = C(3,-gnssRPY(2,3))*C(2,-gnssRPY(2,2))*C(1,-gnssRPY(2,1));
Cpe_1 = Cne_1*Cbn_1;
Omega_iep = inv(Cpe_1)*Omega_iee*Cpe_1;
w_ieP(1,1) = Omega_iep(3,2);
w_ieP(2,1) = Omega_iep(1,3);
w_ieP(3,1) = Omega_iep(2,1);
  
%Quaternion t0 im e-system berechnen 
qt0 = compact(quaternion((Cpe_0),'rotmat','frame'));

%Quaternion t1 im e-system berechnen
qt1 = compact(quaternion((Cpe_1),'rotmat','frame'));
q_e(1,:) = qt0;
q_e(2,:) = qt1;


 %% InertialNavi only RK3

% Rk3_temp = [refp(1,:),refv(1,:),q_e(1,:)];
% for i = 2:length(t_imu)
%     delta_t = abs(t_imu(i) - t_imu(i-1));
%     if i == 2
%         x0 = [IMUdata(i-1, 1:6); IMUdata(i, 1:6)];
%     else
%         x0 = [IMUdata(i-1:i, 1:3) IMUdata(i-1:i, 4:6)];
%     end
%     Rk3_temp(i,:) = RungeKutta3(@TimeDerivativePosVelAtt_e, Rk3_temp(i-1,:)', x0', delta_t);
% 
% 
% end
% RK3_Pos = Rk3_temp(:,1:3);
% RK3_Vel = Rk3_temp(:,4:6);
% RK3_q = Rk3_temp(:,7:10);
% 
% for i = 1:length(RK3_q)
%     RK3_rotm{i,1} = rotmat(quaternion(RK3_q(i,:)),'frame');
% end
% 
% 
% 
% figure
% plot3(gnssp_e(:,1),gnssp_e(:,2),gnssp_e(:,3))
% hold on
% plot3(RK3_Pos(:,1),RK3_Pos(:,2),RK3_Pos(:,3))
% plot3(refp(:,1),refp(:,2),refp(:,3))
% xlabel("[m]");ylabel("[m]");zlabel("[m]");
% title('Aufgabe 1 3D Plot for all');
%% prepare for KF

updaterate = 20; % 1 oder 5 oder 20, [s]
t_Update = t_gnss(1):updaterate:t_gnss(end);

bias_ap = ones(3,1)*1e-12; % initialer Bias Accelerometer
bias_wpip = ones(3,1)*1e-12; % initialer Bias Gyro
Rn = diag([2.0 4.0 3.0 2*1e-1 4*1e-1 2*1e-1].^2); % Messrauschen/Unsicherheit Messungen (Pos [m] und Vel [m/s])
Hn = [eye(6), zeros(6,9)]; % Desingmatrix fuer Beobachtungen
xnn(1,:) = ones(15,1)*1e-12; % initialer Zustandsvektor der Fehler [dPos, dVel, dOri, Bias Acc, Bias Gyro]'
Pnn{1} = diag([ones(1,3)*1e-1...            
               ones(1,3)*1e-1...           
               deg2rad(ones(1,3)*1e-2)...  
               ones(1,3)*1e-1...           
               deg2rad(ones(1,3)*1)].^2);  % zugehoerige Kovarianzmatrix
xe(1,:) = gnssp_e(1,:);
ve(1,:) = gnssv_e(1,:);
Cpe = Cpe_0;

for i = 1:length(t_imu)-1    

   delta_t = abs(t_imu(i) - t_imu(i+1));
   idx_GNSS = find(t_imu(i) == t_gnss);
    
   ae = Cpe*IMUdata(i,1:3)';
   Ae = ome2Ome(ae);
   tao = 3000; %更换其他数据时更换
   beta = 1/tao;
   Fa = -beta*eye(3);
   Fw = -beta*eye(3);
   F = [zeros(3,3) eye(3) zeros(3,3) zeros(3,3) zeros(3,3)
     -Omega_iee*Omega_iee -2*Omega_iee Ae Cpe_0 zeros(3)
     zeros(3,3) zeros(3,3) -Omega_iee zeros(3,3) -Cpe_0
     zeros(3,3) zeros(3,3) zeros(3,3) Fa zeros(3,3)
     zeros(3,3) zeros(3,3) zeros(3,3) zeros(3,3) Fw];
   varianz_ap = 1e-4^2;    % Varianz Accelerometer Bias Rauschen [m^2/s^4]
   varianz_wpip = 1e-5^2;  % Varianz Gyro Bias Rauchen [rad^2/s^2]
   Ga = sqrt(varianz_ap*beta)*eye(3); 
   Gw = sqrt(varianz_wpip*beta)*eye(3); 
   G = [zeros(9,6);
     Ga zeros(3,3);
    zeros(3,3) Gw];
   W = 1;
   A = [-F G*W*G'; zeros(size(F)) F']*mean(diff(t_gnss));
   n = length(A)/2;
   B = expm(A);
   Phi = B(n+1:2*n,n+1:2*n)';  % Zustandsuebergangsmatrix
   Q = Phi*B(1:n,n+1:2*n);     % Matrix des Prozessrauschens
      
    if ~isempty(idx_GNSS) && any(t_gnss(idx_GNSS) == t_Update) % ja, Beobachtungen vorhanden
        
        zn = [gnssp_e(idx_GNSS,:)-xe(i,:) ...  % Positionsdifferenz
              gnssv_e(idx_GNSS,:)-ve(i,:)]';   % Geschwindigkeitsdifferenz

        [x, P] = KalmanFilter(xnn(i,:)', Pnn{i} ,Phi, zn, Hn, Q, Rn);
        xnn(i+1,:) = x;
        Pnn{i+1,:} = P;
        sigma(i+1,:) = sqrt(diag(P));

        % geschaetzte Fehler als Korrekturen anbringen (vor Integration
        % (VO06, F7))
        xe(i,:) = xe(i,:) + xnn(i+1,1:3);   % Position
        ve(i,:) = ve(i,:) + xnn(i+1,4:6);   % Geschwindigkeit
            psi = xnn(7:9);
            Psi = [0 -psi(3) psi(2); psi(3) 0 -psi(1); -psi(2) psi(1) 0];
        Cpe = (eye(3) - Psi)*Cpe;     % Orientierung
        q_e(i,:) = compact(quaternion((Cpe),'rotmat','frame'));   % Orientierung
        
        bias_ap = bias_ap+xnn(i+1,10:12)';      % akkumulierter Bias
        bias_wpip = bias_wpip+xnn(i+1,13:end)'; % akkumulierter Bias

        % Fehler zuruecksetzen (Zuschlaege ab jetzt wieder Null, da
        % Korrektur angebracht ist)
        xnn(i+1,1:6) = xnn(1,1:6); 
        
    else
        % Bestimme Fehler fuer den Fall, dass keine Beob. vorliegen
        % (Praediktion Fehlervektor - als Initialwerte fuer nächste Korrektur
        % mit KalmanFilter)
        zn = [];
        [x, P] = KalmanFilter(xnn(i,:)', Pnn{i} ,Phi, zn, Hn, Q, Rn);
        xnn(i+1,:) = x;
        Pnn{i+1,:} = P;
        sigma(i+1,:) = sqrt(diag(P));
    end

    % (korrigierte) Startwerte RungeKutta
    Rk3_KF(i,:) = [xe(i,:),ve(i,:),q_e(i,:)];
    % Anbringen des geschaetzten Bias' an die Messungen
    IMUdataKF(i,1:3) = IMUdata(i,1:3)+bias_ap';
    IMUdataKF(i,4:6) = IMUdata(i,4:6)+bias_wpip';
    if i == 1
        x0 = [IMUdataKF(i, 1:6); IMUdataKF(i, 1:6)];
    else
        x0 = [IMUdataKF(i-1:i, 1:3) IMUdataKF(i-1:i, 4:6)];
    end
    % Integration e-System (RungeKutta 3.Ordnung)
    Rk3_KF(i+1,:) = RungeKutta3(@TimeDerivativePosVelAtt_e_g, Rk3_KF(i,:)', x0', delta_t);

    % Extrahiere Werte aus Integration
    xe(i+1,:) = Rk3_KF(i,1:3);    % Position ECEF [m m m]
    ve(i+1,:) = Rk3_KF(i,4:6);    % Geschwindigkeit ECEF [m/s m/s m/s]
    q_e(i+1,:) = Rk3_KF(i,7:10);         % Quaternion p2e-System
    Cpe  = rotmat(quaternion(q_e(i+1,:)),'frame');% DCM p2e-System
end

figure
plot3(gnssp_e(:,1),gnssp_e(:,2),gnssp_e(:,3))
hold on
plot3( xe(:,1), xe(:,2), xe(:,3))
plot3(refp(:,1),refp(:,2),refp(:,3))
xlabel("[m]");ylabel("[m]");zlabel("[m]");
legend('gnss','RK3','REF')
title('3D Plot for all');



figure
plot(gnssp_e(:,2),gnssp_e(:,3))
hold on 
plot(xe(:,2),xe(:,3))
plot(refp(:,2),refp(:,3))
xlabel("[m]");ylabel("[m]");zlabel("[m]");
legend('gnss','RK3','REF')
title('2D Plot for all');



figure
subplot(3, 1, 1);
plot(t_gnss,gnssv_n(:,1));
hold on
plot(t_imu,ve(:,3))
plot(t_imu,refv_n(:,1))
xlabel("[s]");ylabel("[m/s]");
legend('gnss','RK3','REF')
title('Geschwindigkeit in N');

subplot(3, 1, 2);
plot(t_gnss,gnssv_n(:,2));
hold on
plot(t_imu,ve(:,2))
plot(t_imu,refv_n(:,2))
ylabel("[m/s]");
legend('gnss','RK3','REF')
title('Geschwindigkeit in E');

subplot(3, 1, 3);
plot(t_gnss,gnssv_n(:,3));
hold on
plot(t_imu,ve(:,1))
plot(t_imu,refv_n(:,3))
ylabel("[m/s]");
legend('gnss','RK3','REF')
title('Geschwindigkeit in D');



