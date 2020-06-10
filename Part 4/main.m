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
spline_order            = 4;
spline_continuity       = 1;

% Model Quality Settings
max_spline_order        = 4;
max_continuity          = max_spline_order - 1;
max_simplices_xy        = 8;

% Plotting Settings
plot_spline     = 1;
plot_validation = 1;                 % This takes a while if turned on!

% Triangulation Settings
num_triangles_x = 4;
num_triangles_y = 4;

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

%% Run Simplex Spline
% Fixed polynomial degree, continuity and number of simplices
[global_B_id, global_B_val, global_idx_val Y_hat_spline, c_spline, VAR,...
    T, x, y, vertices, RMSE] = simplex_continuity1(spline_order,...
    spline_continuity,num_triangles_x, num_triangles_y, X_id, Y_id,...
    X_val, Y_val, plot_spline, save);

% Simplex Spline Plots & Validation of Model
spline_plot(spline_order, spline_continuity, X_id, Y_id, X_val, Y_val,...
    Y_hat_spline, global_idx_val, T, x, y, vertices, RMSE,...
    num_triangles_x, num_triangles_y, plot_spline, save);
validation_spline(spline_order, spline_continuity, num_triangles_x,...
    num_triangles_y, X_id, Y_id, X_val, Y_val, max_spline_order, max_continuity, ...
    max_simplices_xy, Y_hat_spline, global_B_val, global_idx_val, c_spline,...
    VAR, plot_spline, plot_validation, save);
