function qdot = TimeDerivativeQuaterionen(omega_epp,q_pe)
w1 = omega_epp(1);
w2 = omega_epp(2);
w3 = omega_epp(3);

A = [0 w1 w2 w3
    -w1 0 w3 -w2
    -w2 -w3 0 w1
    -w3 w2 -w1 0];
qdot = (1/2*A*q_pe')';
end