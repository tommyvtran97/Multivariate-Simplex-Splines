% SIMPLEX_POLYNOMIAL creates the structure for the simplex polynomial and
% uses an ordinary least square regression to calculate the estimated
% parameters for the simplex polynomial

function [TRI, PHI, Bx_val, c_hat, X_val, Y_val, Yb_hat_val, residual, RMSE] = ...
    simplex_polynomial(X_id, Y_id, X_val, Y_val, simplex_order, plot_result, save)
    
    % Initalize measurement vectors
    alpha_m = X_id(1,:);
    beta_m  = X_id(2,:);
    X_id    = X_id';
    X_val   = X_val';
    
    % Define simplex vertices
    PHI = [0.8 -0.2 ; -0.2 -0.2; -0.2 0.2];
    
    % Define simplex
    TRI = delaunayn(PHI);
    
    % Calculate points inside simplex and barocentric coordinates
    [IMap_id, BaryC_id]    = tsearchn(PHI,TRI,X_id);
    [IMap_val, BaryC_val]  = tsearchn(PHI,TRI,X_val);
  
    IMap_id  = find(~isnan(IMap_id));   IMap_val  = find(~isnan(IMap_val)); 
    BaryC_id = BaryC_id(IMap_id, :);    BaryC_val = BaryC_val(IMap_val, :);
    Y_id     = Y_id(:, IMap_id)';       Y_val     = Y_val(:, IMap_val)';
    X_id     = X_id(IMap_id, :);        X_val     = X_val(IMap_val, :);
    
    % Create B regression matrix
    [B_sorted] = sorted_bcoefficient(simplex_order);
    
    Bx_id      = (zeros(size(BaryC_id, 1), size(B_sorted, 1)));
    Bx_val     = (zeros(size(BaryC_val, 1), size(B_sorted, 1)));
    
    for i=1:1:size(Bx_id, 1)
        for j=1:1:size(Bx_id,2)
            % Define the multinomial
            multi_nom = factorial(simplex_order) / (factorial(B_sorted(j,1)) ...
                        * factorial(B_sorted(j,2)) * factorial(B_sorted(j,3)));
                    
            Bx_id(i,j) = multi_nom * BaryC_id(i,1)^(B_sorted(j,1)) * ...
                BaryC_id(i,2)^(B_sorted(j,2)) * BaryC_id(i,3)^(B_sorted(j,3));
        end
    end
    
    for i=1:1:size(Bx_val, 1)
        for j=1:1:size(Bx_val,2)
            % Define the multinomial
            multi_nom = factorial(simplex_order) / (factorial(B_sorted(j,1)) ...
                        * factorial(B_sorted(j,2)) * factorial(B_sorted(j,3)));
                    
            Bx_val(i,j) = multi_nom * BaryC_val(i,1)^(B_sorted(j,1)) * ...
                BaryC_val(i,2)^(B_sorted(j,2)) * BaryC_val(i,3)^(B_sorted(j,3));
        end
    end
    
    % Calculate the B-coefficients 
    c_hat = pinv(Bx_id' * Bx_id) * Bx_id' * Y_id;
    
    % Calculate estimation from the validation dataset
    Yb_hat_val  = (Bx_val * c_hat)';    % Tranpose to correct size
    
    X_val = X_val';
    Y_val = Y_val';
    
    % Calculate the root mean square (RMS)
    residual = (Y_val - Yb_hat_val);
    RMSE = (rms(residual) / (max(Y_val) - min(Y_val)))*100;

end 