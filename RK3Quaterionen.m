function [rotm,quat] = RK3Quaterionen(func,omega,q0,h)

k1 = func(omega(1,:),q0);
k2 = func(omega(2,:),q0+k1*h/2);
k3 = func(omega(3,:),q0-k1*h+k2*2*h);
q = normalize(quaternion(q0 + h/6*(k1+4*k2+k3)));
quat = compact(q);
rotm = rotmat(q,'frame');
end