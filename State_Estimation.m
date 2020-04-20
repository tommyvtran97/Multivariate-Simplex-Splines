%% State and Parameter Estimation with F16 Flight Data
clc;
clear;

dataname = 'Dataset/F16traindata_CMabV_2020';
load(dataname, 'Cm', 'Z_k', 'U_k');

save = 1;
doIEKF = 0;

perc_val = 0.5;

% Transpose the matrices
Cm = Cm'; Z_k = Z_k'; U_k = U_k';

% measurements Z_k = Z(t) + v(t)
alpha_m = Z_k(1,:); % measured angle of attack
beta_m = Z_k(2,:);  % measured angle of sideslip
Vtot = Z_k(3,:);    % measured velocity

% input to Kalman filter
Au = U_k(1,:); % perfect accelerometer du/dt data
Av = U_k(2,:); % perfect accelerometer dv/dt data
Aw = U_k(3,:); % perfect accelerometer dw/dt data

% Standard deviation for the system noise and measurement noise statistics
stdw    = [1e-3, 1e-3, 1e-3, 0];     % u, v, w and C
stdv    = [0.035, 0.013, 0.110];     % alpha, beta and V

%% Set Initial Values for States and Statistics
Ex_0    = [150; 0.0; 0.0; 0.0];  % initial estimate of optimal value of x_k_1k_1
n       = length(stdw);         % number of states
nm      = size(Z_k,1);          % number of measurements
m       = size(U_k,2);          % number of inputs

B       = eye(n);               % input matrix
G       = eye(length(stdw));    % noise input matrix

%% Initial Estimate for Covariance Matrix
stdx_0  = [0.1, 0.1, 0.1, 0.1];     % Convergence KF depends on this estimate!
P_0     = diag(stdx_0.^2);          % Create diagonal covariance matrix

%% System Noise and Measurement Noise statistics
Q = diag(stdw.^2);
R = diag(stdv.^2);
%w_k = stdw * randn(n, N) + Ew;  %TODO: Is this required?
%v_k = stdv * randn(n, N) + Ev;  %TODO: Is this required?

x_k_1k_1 = Ex_0;            % x(0|0)=E{x_0}
P_k_1k_1 = P_0;             % P(0|0)=P(0)

%% Check for Observability prove that the Kalman Filter converges
%observability();

%% Run Iterated Extended Kalman filter
[z_pred, XX_k1k1, dt] = IEKF_function(U_k, Z_k, stdw, stdv, doIEKF);
%IEKF_plot(alpha_m, beta_m, Vtot, z_pred, U_k, XX_k1k1, save)

%% Split Data into Idenfication and Validation Set
[X_id, X_val, Y_id, Y_val] = split_data(z_pred, Cm, perc_val);

%% Run Ordinary Least Square Estimator
expo = create_polynomial(2, 3);
[Y_hat_id, Y_hat_val, Y_hat_full, theta_hat_id, theta_hat_val, theta_hat_full, A_matrix_id, A_matrix_val, A_matrix_full] = OLS_function(X_id, Y_id', X_val, Y_val',z_pred, Cm', expo);
OLS_plot(X_val, Y_val, Y_hat_val);












