%% State and Parameter Estimation with F16 Flight Data
clc;
clear;

dataname = '../Dataset/F16traindata_CMabV_2020';
load(dataname, 'Cm', 'Z_k', 'U_k');

% Transpose the matrices
Cm = Cm'; Z_k = Z_k'; U_k = U_k';

% Settings
save = 0;
IEKF = 1;

% Multivariate Spline Setting
polynomial_order        = 4;

% Model Quality Settings
max_polynomial_order    = 15;

% Plotting Settings
plot_kalman     = 1;
plot_OLS        = 1;
plot_validation = 1;                 % This takes a while if turned on!

% Standard deviation for the system noise and measurement noise statistics
stdw    = [1e-3, 1e-3, 1e-3, 0];     % u, v, w and C
stdv    = [0.035, 0.013, 0.110];     % alpha, beta and V

%% Check for Observability prove that the Kalman Filter converges
observability();

%% Run Iterated Extended Kalman filter
[z_pred,z_pred_corr, XX_k1k1, IEKFitcount, N_states] = ...
    IEKF_function(U_k, Z_k, stdw, stdv, IEKF);

% IEFK Plots
IEKF_plot(Z_k, z_pred, z_pred_corr, XX_k1k1,...
    IEKFitcount, plot_kalman, save)

%% Split Data into Idenfication and Validation Set
[X_id, X_val, Y_id, Y_val] = split_data(z_pred_corr, Cm);

%% Run Ordinary Least Square Estimator
[Y_hat_val, theta_hat, Ax_val, RMSE] = OLS_function(polynomial_order, X_id, ...
    Y_id, X_val, Y_val);

% OLS Plots & Validation of Model
OLS_plot(polynomial_order, X_val, Y_val, Y_hat_val, RMSE, plot_OLS, save)
validation_polynomial(X_id, Y_id, X_val, Y_val, ...
    polynomial_order, max_polynomial_order, plot_OLS, plot_validation, save)
