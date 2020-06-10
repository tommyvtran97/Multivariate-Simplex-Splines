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
simplex_order           = 10;

% Model Quality Settings
max_simplex_order       = 15;

% Plotting Settings
plot_simplex    = 1;
plot_validation = 1;                 % This takes a while if turned on!

% Standard deviation for the system noise and measurement noise statistics
stdw    = [1e-3, 1e-3, 1e-3, 0];     % u, v, w and C
stdv    = [0.035, 0.013, 0.110];     % alpha, beta and V

%% Check for Observability prove that the Kalman Filter converges
observability();

%% Run Iterated Extended Kalman filter
[z_pred,z_pred_corr, XX_k1k1, IEKFitcount, N_states] = ...
    IEKF_function(U_k, Z_k, stdw, stdv, IEKF);

%% Split Data into Idenfication and Validation Set
[X_id, X_val, Y_id, Y_val] = split_data(z_pred_corr, Cm);

%% Run Simplex Polynomial
[TRI, PHI, Bx_val, c_hat, X_val_simplex, Y_val_simplex, Yb_hat_val, residual,...
    RMSE] = simplex_polynomial(X_id, Y_id, X_val, Y_val, simplex_order,...
    plot_simplex, save);

% Simplex Plots & Validation of Model
simplex_plot(TRI, PHI, X_id, Y_id, X_val_simplex, Y_val_simplex, Yb_hat_val,...
    simplex_order, RMSE, plot_simplex, save);
validation_simplex(X_id, X_val_simplex, Y_id, Y_val_simplex, c_hat,...
    simplex_order, max_simplex_order, plot_simplex, plot_validation, save)