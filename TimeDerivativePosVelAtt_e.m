function [y_dot] = TimeDerivativePosVelAtt_e(y,x)
   w_E = [0;0;7.292115e-5];
   Omega_iee = [0 -w_E(3) 0
             w_E(3) 0 0
             0 0 0];
 
  x_current = x(:,1);
  x_previous = x(:,2);
  % Data extraction
  pos_e = y(1:3);
  v_e = y(4:6);
  
  C_p2e_quat = y(7:10);
  C_p2e_rotm = rotmat(quaternion(C_p2e_quat'),'frame');
%这里修改
  a_ip_p = (x_current(1:3) + x_previous(1:3)) / 2;
  omega_ip_p = (x_current(4:6) + x_previous(4:6)) / 2;

  % Calculations
  omega_ep_p = omega_ip_p - C_p2e_rotm' * w_E;
  A = [      0.0     , omega_ep_p(1), omega_ep_p(2), omega_ep_p(3);...
        -omega_ep_p(1),       0.0     ,  omega_ep_p(3), -omega_ep_p(2);...
        -omega_ep_p(2), -omega_ep_p(3),       0.0     ,  omega_ep_p(1);...
        -omega_ep_p(3),  omega_ep_p(2), -omega_ep_p(1),       0.0     ];

  y_dot = [v_e
           C_p2e_rotm * a_ip_p - 2 * Omega_iee * v_e -...
           Omega_iee * Omega_iee * pos_e
           0.5 * A * C_p2e_quat];
end