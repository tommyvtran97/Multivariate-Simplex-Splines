%% State and Parameter Estimation with F16 Flight Data
clc;
clear;

dataname = 'Dataset/F16traindata_CMabV_2020';
load(dataname, 'Cm', 'Z_k', 'U_k');

% Settings
save = 1;
doIEKF = 1;

% Multivariate Spline Setting
polynomial_order = 10;

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
Ex_0    = [Z_k(3,1); 0.5; 0.5; 0.5];  % initial estimate of optimal value of x_k_1k_1
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
observability();

%% Run Iterated Extended Kalman filter
[z_pred, XX_k1k1, dt, IEKFitcount] = IEKF_function(U_k, Z_k, stdw, stdv, doIEKF);
IEKF_plot(alpha_m, beta_m, Vtot, Cm, z_pred, U_k, XX_k1k1, IEKFitcount, save)

%% Split Data into Idenfication and Validation Set
[X_id, X_val, Y_id, Y_val] = split_data(z_pred, Cm);

%% Run Ordinary Least Square Estimator
MSE_x = [];
MSE_y = [];

% Calculate the Mean Squared Error for different polynomial order
for i=1:1:polynomial_order
    expo = create_polynomial(n, i);
    [Y_hat_id, Y_hat_val, theta_hat, A_matrix_val] = OLS_function(X_id, Y_id, X_val, expo);

    residual  = (Y_val - Y_hat_val).^2;
    MSE = sum(residual)/size(residual, 2);

    MSE_x = [MSE_x, i];
    MSE_y = [MSE_y, MSE];
end

[VAR, COR] = OLS_plot(X_val, Y_val, Y_hat_val, MSE_x, MSE_y, A_matrix_val, theta_hat, save); 












