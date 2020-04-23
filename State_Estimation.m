%% State and Parameter Estimation with F16 Flight Data
clc;
clear;

dataname = 'Dataset/F16traindata_CMabV_2020';
load(dataname, 'Cm', 'Z_k', 'U_k');

% Transpose the matrices
Cm = Cm'; Z_k = Z_k'; U_k = U_k';

% Settings
save = 0;
doIEKF = 1;

% Multivariate Spline Setting
max_polynomial_order    = 20;
polynomial_order        = 10;

% Plotting Settings
plot_kalman = 0;
plot_OLS    = 0;

% Measurements Z_k = Z(t) + v(t)
alpha_m = Z_k(1,:); % measured angle of attack
beta_m = Z_k(2,:);  % measured angle of sideslip
Vtot = Z_k(3,:);    % measured velocity

% Standard deviation for the system noise and measurement noise statistics
stdw    = [1e-3, 1e-3, 1e-3, 0];     % u, v, w and C
stdv    = [0.035, 0.013, 0.110];     % alpha, beta and V

%% Check for Observability prove that the Kalman Filter converges
observability();

%% Run Iterated Extended Kalman filter
[z_pred,z_pred_correct, XX_k1k1, IEKFitcount, N_states] = IEKF_function(U_k, Z_k, stdw, stdv, doIEKF);

% IEFK Plots
IEKF_plot(alpha_m, beta_m, Vtot, Cm, z_pred, z_pred_correct, U_k, XX_k1k1, IEKFitcount, save, plot_kalman)

%% Split Data into Idenfication and Validation Set
[X_id, X_val, Y_id, Y_val] = split_data(z_pred_correct, Cm);

%% Run Ordinary Least Square Estimator
expo = create_polynomial(polynomial_order);
[Y_hat_id, Y_hat_val, theta_hat, Ax_val] = OLS_function(X_id, Y_id, X_val, expo);

% OLS Plots & Validation of Model
OLS_plot(X_val, Y_val, Y_hat_val, save, plot_OLS)
validation(X_id, X_val, Y_id, Y_val, Y_hat_val, polynomial_order, max_polynomial_order, save, plot_OLS)
