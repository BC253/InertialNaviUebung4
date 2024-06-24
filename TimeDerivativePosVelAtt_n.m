function [y_dot] = TimeDerivativePosVelAtt_n(y,x)
   w_E = [0;0;7.292115e-5];
   Omega_iee = [0 -w_E(3) 0
             w_E(3) 0 0
             0 0 0];
  x_current = x(:,1);
  x_previous = x(:,2);

  a_ip_p = (x_current(1:3) + x_previous(1:3)) / 2;
  omega_ipp = (x_current(4:6) + x_previous(4:6)) / 2;
  % Data extraction
  pos_e = lla2ecef([rad2deg(y(1:2)') y(3)]);
  v_n = y(4:6);
  phi = y(1);
  lambda = y(2);
  h = y(3);
  % 2.5
  M = rcurve('meridian',referenceEllipsoid('WGS84'),rad2deg(phi));
  N = rcurve('transverse',referenceEllipsoid('WGS84'),rad2deg(phi));
  phi_dot = v_n(1)/(M+h);
  lambda_dot = v_n(2)/((N+h)*cos(phi));
  h_dot = -v_n(3);
  Cne = C(3,-lambda)*C(2,phi+pi/2);
  % 2.10
  omega_ien = Cne'*w_E;
  omega_enn = [cos(phi)*lambda_dot
               -phi_dot
               -sin(phi)*lambda_dot];
  omega_inn = omega_ien+omega_enn; %Inertialsensorik V06S7
  C_p2n_quat = y(7:10);
  C_p2n_rotm = rotmat(quaternion(C_p2n_quat'),'frame');
  omega_inp = C_p2n_rotm'*omega_inn;
  omega_np_p = omega_ipp - omega_inp;
  
  Omega_ie_n = ome2Ome(omega_ien);
  Omega_en_n = ome2Ome(omega_enn);
  
  g_n = Somigliana(rad2deg(phi),h);
  g_n = [0;0;g_n];

  % A matrix Calculations
  A = [      0.0     , omega_np_p(1), omega_np_p(2), omega_np_p(3);...
        -omega_np_p(1),       0.0     ,  omega_np_p(3), -omega_np_p(2);...
        -omega_np_p(2), -omega_np_p(3),       0.0     ,  omega_np_p(1);...
        -omega_np_p(3),  omega_np_p(2), -omega_np_p(1),       0.0     ];

  y_dot = [phi_dot
           lambda_dot
           h_dot
           C_p2n_rotm * a_ip_p - (2 * Omega_ie_n + Omega_en_n) * v_n -...
           Cne'*Omega_iee*Omega_iee*pos_e'+g_n; %2.11
           0.5 * A * C_p2n_quat];
end