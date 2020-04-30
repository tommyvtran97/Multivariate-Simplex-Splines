%% State and Parameter Estimation with F16 Flight Data
clc;
clear;

dataname = 'Dataset/F16traindata_CMabV_2020';
load(dataname, 'Cm', 'Z_k', 'U_k');

% Transpose the matrices
Cm = Cm'; Z_k = Z_k'; U_k = U_k';

% Settings
save = 0;
IEKF = 1;

% Multivariate Spline Setting
max_polynomial_order    = 15;
polynomial_order        = 10;
max_simplex_order       = 15;
simplex_order           = 10;
max_spline_order        = 3;
max_continuity          = max_spline_order - 1;
spline_order            = 3;
spline_continuity       = 1;

% Plotting Settings
plot_kalman     = 0;
plot_OLS        = 0;
plot_simplex    = 0;
plot_spline     = 1;
plot_validation = 0;

% Triangulation Settings
num_triangles_x = 1;
num_triangles_y = 1;

    
% Standard deviation for the system noise and measurement noise statistics
stdw    = [1e-3, 1e-3, 1e-3, 0];     % u, v, w and C
stdv    = [0.035, 0.013, 0.110];     % alpha, beta and V

%% Check for Observability prove that the Kalman Filter converges
%observability();

%% Run Iterated Extended Kalman filter
[z_pred,z_pred_corr, XX_k1k1, IEKFitcount, N_states] = ...
    IEKF_function(U_k, Z_k, stdw, stdv, IEKF);

% IEFK Plots
IEKF_plot(Z_k, Cm, z_pred, z_pred_corr, U_k,XX_k1k1,...
    IEKFitcount, plot_kalman, save)

%% Split Data into Idenfication and Validation Set
[X_id, X_val, Y_id, Y_val] = split_data(z_pred_corr, Cm);

%% Run Ordinary Least Square Estimator
expo = create_polynomial(polynomial_order);
[Y_hat_val, theta_hat, Ax_val] = OLS_function(X_id, Y_id, X_val, expo);

% OLS Plots & Validation of Model
OLS_plot(X_val, Y_val, Y_hat_val, plot_OLS, save)
validation_polynomial(X_id, Y_id, X_val, Y_val, ...
    polynomial_order, max_polynomial_order, plot_validation, save)

%% Run Simplex Polynomial
[TRI, PHI, Bx_val, c_hat, X_val_simplex, Y_val_simplex, Yb_hat_val, residual, RMSE] = ...
    simplex_polynomial(X_id, Y_id, X_val, Y_val, simplex_order,...
    plot_simplex, save);

% Simplex Plots & Validation of Model
simplex_plot(TRI, PHI, X_id, Y_id, X_val_simplex, Y_val_simplex, Yb_hat_val,...
    simplex_order, plot_simplex, save);
    
validation_simplex(X_id, X_val_simplex, Y_id, Y_val_simplex, c_hat,...
    simplex_order, max_simplex_order, plot_validation, save)

%% Run Simplex Spline

% Fixed continuity and varying polynomial degree model performance
[global_B_id, global_B_val, global_idx_val Y_hat_spline, c_spline, VAR, RMSE_x, RMSE_y, H, Tri, T] = simplex_continuity1(spline_order, spline_continuity,...
    num_triangles_x, num_triangles_y, X_id, Y_id, X_val, Y_val, max_spline_order, plot_spline);

% Fixed polynomial degree and varying continuity model performance
[global_B_id, global_B_val, global_idx_val Y_hat_spline, c_spline, VAR, RMSE_x_cont, RMSE_y_cont, H, Tri, T] = simplex_continuity2(spline_order, max_continuity,...
    num_triangles_x, num_triangles_y, X_id, Y_id, X_val, Y_val, max_spline_order, plot_spline);

% Fixed continuity and polynomial degree
[global_B_id, global_B_val, global_idx_val Y_hat_spline, c_spline, VAR, H, Tri, T] = simplex_continuity3(spline_order, spline_continuity,...
    num_triangles_x, num_triangles_y, X_id, Y_id, X_val, Y_val, max_spline_order, plot_spline);

% Simplex Spline Plots & Validation of Model
spline_plot(spline_order, X_id, Y_id, X_val, Y_val,...
    Y_hat_spline, plot_spline, save);
validation_spline(spline_order, Y_hat_spline, ...
    Y_val,global_B_val, global_idx_val, c_spline, VAR, RMSE_x, RMSE_y, RMSE_x_cont, RMSE_y_cont, plot_validation, save);
