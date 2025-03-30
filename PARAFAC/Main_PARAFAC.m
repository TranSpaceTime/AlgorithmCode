clf;clc; clear; close all;

%% 参数设置
M = 6;           % x轴阵元个数
N = 5;           % y轴阵元个数
K = 512;        % 快拍数
fc = 100e6;     % 载波频率
fs = 300e6;     % 采样频率
Pn = 1;          % 噪声功率

fines = [20, 40]; % 方位角（度）
thetas = [5, 60];    % 俯仰角（度）
signal_SNR = [30, 30]; % 信噪比（dB）
signal_f = [15e6, 30e6]; % 信号频率

m = (0:M-1)';    % x轴阵元索引
n = (0:N-1)';    % y轴阵元索引
c = 3e8;         % 光速
lambda = c / fc;  % 波长
dx = lambda / 2;  % x轴阵元间距
dy = lambda / 2;  % y轴阵元间距
num_signals = length(fines); % 信号数目

%% 生成信号
t = (0:K-1)/fs; % 时间向量
S = zeros(num_signals, K);
for k = 1:num_signals
    A_k = sqrt(Pn)*10^(signal_SNR(k)/20); % 信号幅度
    S(k,:) = A_k * exp(1j*2*pi*signal_f(k)*t); % 载波调制
end

%% 构造阵列流型
A = zeros(M, num_signals);
B = zeros(N, num_signals);
for k = 1:num_signals
    phi = deg2rad(fines(k));
    theta = deg2rad(thetas(k));
    
    % 空间频率计算
    u = (dx/lambda)*sin(theta)*cos(phi);
    v = (dy/lambda)*sin(theta)*sin(phi);
    
    A(:,k) = exp(-1j*2*pi*m*u); % x方向导向矢量
    B(:,k) = exp(-1j*2*pi*n*v); % y方向导向矢量
end

%% 构造接收张量
X = tensor(zeros(M,N,K));
noise = (randn(M,N,K) + 1j*randn(M,N,K)) * sqrt(Pn/2);

for k = 1:num_signals
    component = ktensor(1, A(:,k), B(:,k), S(k,:).');
    X = X + tensor(component);
end
X = X + tensor(noise); 
X_normalized = X / norm(X); % 张量归一化

%% 稳健PARAFAC分解
R = num_signals; % 分解秩
options = struct;
options.init = 'nvecs';  % 使用nvecs初始化
options.printitn = 1;    % 显示迭代过程
options.tol = 1e-4;      % 设置收敛容差
options.maxiters = 150;  % 增加最大迭代次数

[Factors, ~] = cp_als(X_normalized, R, options);

A_est = Factors{1}; 
B_est = Factors{2};

%% 参数估计
estimated_angles = zeros(num_signals, 2);
for d = 1:R
    % x方向估计
    a = A_est(:,d);
    phase_x = unwrap(angle(a));
    u_est = -(phase_x(2:end) - phase_x(1:end-1))/(2*pi*(dx/lambda));
    u_est_avg = mean(u_est);
    
    % y方向估计
    b = B_est(:,d);
    phase_y = unwrap(angle(b));
    v_est = -(phase_y(2:end) - phase_y(1:end-1))/(2*pi*(dy/lambda));
    v_est_avg = mean(v_est);
    
    % 角度解算
    phi_est = -atan2d(v_est_avg , u_est_avg);
    theta_est = asind(sqrt(u_est_avg.^2 + v_est_avg.^2));
    
    estimated_angles(d,:) = [phi_est, theta_est];
end

%% 二维可视化
figure(10);

% 真实位置三维散点图
scatter(fines, thetas, 'r', 'filled');hold on;
% 估计位置三维散点图
scatter(estimated_angles(:,1), estimated_angles(:,2), 'b', 'd');
xlabel('方位角 (度)'); ylabel('俯仰角 (度)'); zlabel('幅度');
title('信号分布');
axis([0,90,0,90]);
grid on; 

%% 性能评估
RMSE_phi = sqrt(mean((fines - estimated_angles(:,1)').^2));
RMSE_theta = sqrt(mean((thetas - estimated_angles(:,2)').^2));
disp(['方位角RMSE: ', num2str(RMSE_phi), ' 度']);
disp(['俯仰角RMSE: ', num2str(RMSE_theta), ' 度']);

