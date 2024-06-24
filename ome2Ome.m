function Omega = ome2Ome(omega)
ome1 = omega(1);ome2 = omega(2);ome3 = omega(3);
Omega = [ 0  -ome3  ome2
         ome3  0   -ome1
        -ome2 ome1   0];
end