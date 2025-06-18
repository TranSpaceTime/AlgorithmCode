%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 基于四阶累积量
% MUSIC-like算法，三种常见方法的测试
% 均匀阵和稀疏阵 均可使用 改变"idx"即可
% 解释：https://blog.csdn.net/qq_63978316/article/details/148486692
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; clf;

%% 参数设定
lambda = 1500;                   % 波长
d = lambda/2;                    % 阵元间距
theta = [-50,-20,0,60];          % 入射角度（单位：度）   -20,30 ||  -50,-20,0,60
theta_rad = theta * pi/180;      % 转为弧度
idx = [0,1,2]';             	   % 阵元索引，阵元
M = length(idx);                 % 阵元数
K = length(theta);               % 信号源数
T = 50000;                       % 快拍数
SNR = 50;                        % 信噪比

%% 信号模型
% S = ( (2*randi([0 1],K,T)-1) + 1j*(2*randi([0 1],K,T)-1) ); % K个BPSK信源
S = generateQAM(2, K, T);  % 生成 K 路 16-QAM 信号
A = exp(-1j*pi*idx*sin(theta_rad));                  % 阵列导向向量矩阵（M×K）
X_clean = A * S;                                     % 干净信号 M×T
signal_power = mean(abs(X_clean(:)).^2);             % 信号平均功率
sigma2 = signal_power / (10^(SNR/10));               % 根据SNR计算噪声功率
N = sqrt(sigma2/2) * (randn(M,T) + 1j*randn(M,T));   % AWGN复噪声
X = X_clean + N;                                     % 加噪后的接收信号

%% 构造四阶累积量的各种方法
% 1.利用概念构造严格四阶累积量矩阵(成功估计)
[C4_funO] = BasicDefinitionCumulant4(X);  % 使用整体观测信号计算四阶累积量 (M^2, M^2)
% 2.利用向量的形式构造四阶累积量矩阵(目前 成功 估计) 
[C4_funT] = SingleSampleCumulant4(X);   
% 3.利用所有样本一同构造四阶累积量矩阵(目前 失败 估计) 
% [C4_funX] = WholeSampleCumulant4(X);
%% 特征值分解，获取噪声子空间
[U,D] = eig(C4_funO);                                 % 特征值分解
[D,I] = sort(diag(D));                                % 将特征值排序从小到大
U = fliplr(U(:, I));                                  % 对应特征矢量排序，fliplr 之后，较大特征值对应的特征矢量在前面
Un = U(:, K+1:M^2);                                   % 噪声子空间
%% MUSIC-Like 空间谱扫描
angle_scan = -90:0.1:90;
Pmusic = zeros(1, length(angle_scan));
for i = 1:length(angle_scan)
    a = exp(-1j*pi*idx*sin(angle_scan(i)*pi/180)); % 导向向量
    b = kron(a, conj(a));                          % 扩展向量（M²×1）
    Pmusic(i) = 1 / real(b' * (Un * Un') * b);     % 空间谱
end

%% 归一化 & 绘图
Pmusic_dB = 10 * log10(abs(Pmusic) / max(abs(Pmusic)));
z_min = min(Pmusic_dB(:));    z_max = max(Pmusic_dB(:));
figure(1);
plot(angle_scan, Pmusic_dB, 'LineWidth', 1.5);
xlabel('\theta (°)', 'FontSize', 12);
ylabel('Spatial Spectrum (dB)', 'FontSize', 12);
grid on;
title('MUSIC-like Spatial Spectrum ');
hold on;
for k = 1:K  % 只显示前K个
    plot([theta(k), theta(k)],[z_min, z_max], '-.', 'LineWidth', 0.5);
end

%% 峰值搜索
[pks, locs] = findpeaks(Pmusic_dB);
[~, idx_max] = maxk(pks, K);
Theta_est = sort(angle_scan(locs(idx_max)));
disp('估计角度为：');
disp(Theta_est);

%% ========== 子函数定义 ==========
function [C] = BasicDefinitionCumulant4(X) 
    % X: M × T 的复信号矩阵
    [M, ~] = size(X);
    C = zeros(M^2, M^2);   
    for m = 1:M
        xm = X(m,:);
        for n = 1:M
            xn = conj(X(n,:));
            for k = 1:M
                xk = conj(X(k,:));
                for l = 1:M
                    xl = X(l,:);
                    E1 = mean(xm .* xn .* xk .* xl);
                    E2 = mean(xm .* xn) * mean(xk .* xl); 
                    E3 = mean(xm .* xk) * mean(xn .* xl); 
                    E4 = mean(xm .* xl) * mean(xn .* xk); % 可删去“充分对称分布”-> 任何非共轭的二阶矩都等于 0                   
                    left  = n + (m-1)*M;
                    right = l + (k-1)*M;                    
                    C(left, right) = E1-E2-E3-E4;                          
                end
            end
        end
    end
end

function [C] = SingleSampleCumulant4(X)
    [M, T] = size(X);   Z = zeros(M^2, T);
    % 计算E1
    for t = 1:T
        Z(:, t) = khatrirao(X(:, t), conj(X(:, t)));
    end
    E1 = 1/T*(Z*Z') ;
    % 计算E2
    E2 = mean(Z, 2) * mean(Z, 2)';  %        
    % 计算E3
    R = 1/T*(X * X');
    E3 = kron(R, conj(R)); 
    % 构造累积量
    C = E1 - E2 - E3; % 
end

function [C] = WholeSampleCumulant4(X) 
    [~, T] = size(X);
    E1 = (kron(X,conj(X))*(kron(X,conj(X))'))/T;
    E2 =  kron(X,conj(X))/T*kron(X,conj(X))' /T;
    E3 =  kron(X*X'/T,conj(X*X'/T));
    C = E1 -E2 -E3;  % 四阶累积量计算            
end

function S = generateQAM(L, K, T)

    % L: 每维星座点数为 2L（例如 L=2 → ±1, ±3 → 16-QAM）
    % K: 信源个数
    % T: 快拍数（时间采样点）
    levels = -(2*L-1):2:(2*L-1);  % 构造对称星座点，例如 ±1, ±3, ±5...    
    % 均匀从 levels 中随机采样实部和虚部
    I = datasample(levels, K*T, 'Replace', true);
    Q = datasample(levels, K*T, 'Replace', true);
    % 构造复数星座
    S = reshape(I + 1j*Q, K, T);    
    % 若 L=1 ⇒ 星座为 ±1+j·±1（4-QAM ≈ QPSK）
    % 若 L=2 ⇒ 星座为 ±1,±3+j·±1,±3（16-QAM）
    % 若 L=4 ⇒ 64-QAM   
    
end
