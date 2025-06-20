%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Fourth‑Order Cumulant Based
% MUSIC-like算法，三种常见方法的测试
% MUSIC‑like algorithm – test of three typical implementations
% 均匀阵和稀疏阵 均可使用 改变"idx"即可
% Works for both uniform and sparse arrays; change "idx" as needed
% My interpretation: https://blog.csdn.net/qq_63978316/article/details/148486692
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; clf;

%% Parameter settings
lambda = 1500;                   % Wavelength
d = lambda/2;                    % Array element spacing
theta = [-50,-20,0,60];          % Incident angles (degrees)   -20,30 ||  -50,-20,0,60
theta_rad = theta * pi/180;      % Convert to radians
idx = [0,1,2]';                  % Element indices
M = length(idx);                 % Number of sensors
K = length(theta);               % Number of sources
T = 50000;                       % Snapshots
SNR = 50;                        % Signal‑to‑noise ratio (dB)

%% Signal model
% S = ( (2*randi([0 1],K,T)-1) + 1j*(2*randi([0 1],K,T)-1) ); % K BPSK sources
S = generateQAM(2, K, T);                 % Generate K 16‑QAM sources
A = exp(-1j*pi*idx*sin(theta_rad));       % Steering matrix (M×K)
X_clean = A * S;                          % Noise‑free received data (M×T)
signal_power = mean(abs(X_clean(:)).^2);  % Average signal power
sigma2 = signal_power / (10^(SNR/10));    % Noise power from desired SNR
N = sqrt(sigma2/2) * (randn(M,T) + 1j*randn(M,T));  % Complex AWGN
X = X_clean + N;                          % Noisy data

%% Build fourth‑order cumulant matrices
% 1. Strict fourth‑order cumulant by definition  (works)
[C4_funO] = BasicDefinitionCumulant4(X);      % (M^2 × M^2)
% 2. Vector‑form fourth‑order cumulant (works)
[C4_funT] = SingleSampleCumulant4(X);   
% 3. Cumulant via entire sample (currently fails)
% [C4_funX] = WholeSampleCumulant4(X);

%% Eigen‑decomposition to obtain noise subspace
[U,D] = eig(C4_funO);                % Eigendecomposition
[D,I] = sort(diag(D));               % Sort eigenvalues ascending
U = fliplr(U(:, I));                 % Rearrange eigenvectors (largest first)
Un = U(:, K+1:M^2);                  % Noise subspace

%% MUSIC‑like spatial spectrum search
angle_scan = -90:0.1:90;
Pmusic = zeros(1, length(angle_scan));
for i = 1:length(angle_scan)
    a = exp(-1j*pi*idx*sin(angle_scan(i)*pi/180));  % Steering vector
    b = kron(a, conj(a));                           % Extended vector (M²×1)
    Pmusic(i) = 1 / real(b' * (Un * Un') * b);      % Spatial spectrum
end

%% Normalise & plot
Pmusic_dB = 10 * log10(abs(Pmusic) / max(abs(Pmusic)));
z_min = min(Pmusic_dB(:));    z_max = max(Pmusic_dB(:));
figure(1);
plot(angle_scan, Pmusic_dB, 'LineWidth', 1.5);
xlabel('\theta (°)', 'FontSize', 12);
ylabel('Spatial Spectrum (dB)', 'FontSize', 12);
grid on;
title('MUSIC‑like Spatial Spectrum');
hold on;
for k = 1:K
    plot([theta(k), theta(k)],[z_min, z_max], '-.', 'LineWidth', 0.5);
end

%% Peak search
[pks, locs] = findpeaks(Pmusic_dB);
[~, idx_max] = maxk(pks, K);
Theta_est = sort(angle_scan(locs(idx_max)));
disp('Estimated DOAs (degrees):');
disp(Theta_est);

%% ========== Function definitions ==========
function [C] = BasicDefinitionCumulant4(X)
    % X: complex data matrix (M × T)
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
                    E4 = mean(xm .* xl) * mean(xn .* xk); % Can be omitted for circularly symmetric sources
                    left  = n + (m-1)*M;
                    right = l + (k-1)*M;
                    C(left, right) = E1 - E2 - E3 - E4;
                end
            end
        end
    end
end

function [C] = SingleSampleCumulant4(X)
    [M, T] = size(X);
    Z = zeros(M^2, T);
    % Compute E1
    for t = 1:T
        Z(:, t) = khatrirao(X(:, t), conj(X(:, t)));
    end
    E1 = 1/T * (Z * Z');
    % Compute E2
    E2 = mean(Z, 2) * mean(Z, 2)';
    % Compute E3
    R = 1/T * (X * X');
    E3 = kron(R, conj(R));
    % Fourth‑order cumulant
    C = E1 - E2 - E3;
end

function [C] = WholeSampleCumulant4(X)
    [~, T] = size(X);
    E1 = (kron(X, conj(X)) * (kron(X, conj(X))')) / T;
    E2 = kron(X, conj(X)) / T * kron(X, conj(X))' / T;
    E3 = kron(X * X' / T, conj(X * X' / T));
    C = E1 - E2 - E3;  % Fourth‑order cumulant
end

function S = generateQAM(L, K, T)
    % L: 每维星座点数为 2L（例如 L=2 → ±1, ±3 → 16-QAM）
    % K: 信源个数
    % T: 快拍数（时间采样点）
    levels = -(2*L-1):2:(2*L-1);          % Symmetric constellation levels
    % Uniform random sampling of I/Q components
    I = datasample(levels, K*T, 'Replace', true);
    Q = datasample(levels, K*T, 'Replace', true);
    % Build complex symbols
    S = reshape(I + 1j*Q, K, T);
    % 若 L=1 ⇒ 星座为 ±1+j·±1（4-QAM ≈ QPSK）
    % 若 L=2 ⇒ 星座为 ±1,±3+j·±1,±3（16-QAM）
    % 若 L=4 ⇒ 64-QAM  
end
