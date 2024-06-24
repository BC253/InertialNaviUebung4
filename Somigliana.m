function g_n = Somigliana(phi,h)

%phi 地理维度（deg） h椭球高
% 定义WGS84常数和参数
a = 6378137; % 长半轴，单位为米
f = 1/298.257223563; % 扁率
b = a * (1 - f); % 短半轴
mu = 3.986004418e14; % 地球的引力常数乘以地球质量，单位为m^3/s^2
omega_e = 7.2921150e-5; % 地球自转角速度，单位为rad/s

% 将纬度从度转换为弧度
phi = deg2rad(phi);

% 计算第一偏心率的平方
e_squared = 2*f - f^2;

% 使用公式计算g0
g0 = 9.7803253359 * (1 + 0.001931853 * sin(phi)^2) / sqrt(1 - e_squared * sin(phi)^2);

% 使用公式计算k
k = 1 - (2*h/a)  * (1 + f + b/ mu*((omega_e^2 * a^2) ))  + 3 * (h^2/a^2);

% 计算重力加速度g
g_n = g0 * k;


end


