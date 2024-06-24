function [y_dot] = TimeDerivativePosVelAtt_e(y,x)
   w_E = [0;0;7.292115e-5];
   Omega_iee = [0 -w_E(3) 0
             w_E(3) 0 0
             0 0 0];
 
  x_current = x(:,1);
  x_previous = x(:,2);
  % Data extraction
  pos_e = y(1:3);
  pos_e_lla = ecef2lla(pos_e');
  phi = pos_e_lla(1);
  h = pos_e_lla(3);
  v_e = y(4:6);
  

 
  C_p2e_quat = y(7:10);
  C_p2e_rotm = rotmat(quaternion(C_p2e_quat'),'frame');
%这里修改
  a_ip_p = (x_current(1:3) + x_previous(1:3)) / 2;
  omega_ip_p = (x_current(4:6) + x_previous(4:6)) / 2;
  g_n = Somigliana(rad2deg(phi),h);
  % g_e = Cne*g_n;
  g_e = -g_n*[pos_e(1)/norm(pos_e),pos_e(2)/norm(pos_e),pos_e(3)/norm(pos_e)]';

  % Calculations
  omega_ep_p = omega_ip_p - C_p2e_rotm' * w_E;
  A = [      0.0     , omega_ep_p(1), omega_ep_p(2), omega_ep_p(3);...
        -omega_ep_p(1),       0.0     ,  omega_ep_p(3), -omega_ep_p(2);...
        -omega_ep_p(2), -omega_ep_p(3),       0.0     ,  omega_ep_p(1);...
        -omega_ep_p(3),  omega_ep_p(2), -omega_ep_p(1),       0.0     ];

   part1 = v_e; % 3x1
   part2 = C_p2e_rotm * a_ip_p - 2 * Omega_iee * v_e - Omega_iee * Omega_iee * pos_e + g_e; % 3x1
   part3 = 0.5 * A * C_p2e_quat; % 4x1

   % Concatenate all parts into y_dot
   y_dot = [part1; part2; part3];
end