function [xnn, Pnn] = KalmanFilter(xnn, Pnn,Phi, zn, Hn, Q, Rn)
% KALMAN FILTER
% Input:    Initialer Zustandsvektor                xnn
%           Initiale Kovarianzmatrix                Pnn
%           Zustandsuebrgangsmatrix                 Phi
%           Beobachtungsvektor                      zn
%           Designmatrix der Beoabchtungen          Hn
%           Kovarianzmatrix des Prozessrauschens    Q
%           Varianzmatrix der Messunsicherheiten    R
% Output:   geschaetzter Zustandsvektor             xnn
%           zugehoerige Kovarianzmatrix             Pnn

% Praediktion
xnn_p = Phi*xnn; % praedizierter Zustandsvektor
Pnn_p = Phi*Pnn*Phi' + Q; % zugehoerige Kovarianzmatrix

if (~isempty(zn)) % wenn Beobachtung vorhanden, dann nutze fuer Update
    % Kalman Gain
    Kn = Pnn_p * Hn' * inv(Hn*Pnn_p*Hn'+Rn);

    % Update
    xnn = xnn_p + Kn * (zn-Hn*xnn_p); 
    Pnn = (eye(length(xnn)) - Kn*Hn)*Pnn_p;
else % sonst nur Praediktion des Fehlervektors (damit ordentliche Startwerte fuer KF, wenn Beob. da)
    xnn = xnn_p;
    Pnn = Pnn_p;
end

end