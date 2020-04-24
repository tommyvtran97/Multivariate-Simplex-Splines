function [TRI, PHI, IMap_id, BaryC_id, Bx_val, c_hat, X_val, Y_val, Yb_hat_val, residual, RMSE] = ...
    simplex_polynomial(X_id, Y_id, X_val, Y_val, simplex_order, plot_result, save)
    
    % Initalize measurement vectors
    alpha_m = X_id(1,:);
    beta_m  = X_id(2,:);
    X_id    = X_id';
    X_val   = X_val';
    Y_old = Y_id;
    
    % Define simplex vertices
    PHI = [0.8 -0.2 ; -0.2 -0.2; -0.2 0.2];
    
    % Define simplex
    TRI = delaunayn(PHI);
    
    % Calculate points inside simplex and barocoordinates
    [IMap_id, BaryC_id] = tsearchn(PHI,TRI,X_id);
    [IMap_val, BaryC_val] = tsearchn(PHI,TRI,X_val);
  
    IMap_id = find(~isnan(IMap_id));    IMap_val = find(~isnan(IMap_val)); 
    BaryC_id = BaryC_id(IMap_id, :);    BaryC_val = BaryC_val(IMap_val, :);
    Y_id = Y_id(:, IMap_id)';           Y_val = Y_val(:, IMap_val)';
    X_id = X_id(IMap_id, :);            X_val = X_val(IMap_val, :);
    
    % Create B regression matrix
    [B_sorted] = sorted_bcoefficient(simplex_order);
    
    Bx_id   = (zeros(size(BaryC_id, 1), size(B_sorted, 1)));
    Bx_val  = (zeros(size(BaryC_val, 1), size(B_sorted, 1)));
    
    for i=1:1:size(Bx_id, 1)
        for j=1:1:size(Bx_id,2)
            % Define the multinomial
            multi_nom = factorial(size(B_sorted, 2)) / (factorial(B_sorted(j,1)) ...
                        * factorial(B_sorted(j,2)) * factorial(B_sorted(j,3)));
                    
            Bx_id(i,j) = multi_nom * BaryC_id(i,1)^(B_sorted(j,1)) * BaryC_id(i,2)^(B_sorted(j,2)) * ...
                        BaryC_id(i,3)^(B_sorted(j,3));
        end
    end
    
    for i=1:1:size(Bx_val, 1)
        for j=1:1:size(Bx_val,2)
            % Define the multinomial
            multi_nom = factorial(size(B_sorted, 2)) / (factorial(B_sorted(j,1)) ...
                        * factorial(B_sorted(j,2)) * factorial(B_sorted(j,3)));
                    
            Bx_val(i,j) = multi_nom * BaryC_val(i,1)^(B_sorted(j,1)) * BaryC_val(i,2)^(B_sorted(j,2)) * ...
                        BaryC_val(i,3)^(B_sorted(j,3));
        end
    end
    
    % Calculate the C coefficients 
    c_hat = pinv(Bx_id' * Bx_id) * Bx_id' * Y_id;
    
    % Calculate validation estiation
    Yb_hat_val  = (Bx_val * c_hat)';    % Tranpose to correct size
    
    X_val = X_val';
    Y_val = Y_val';
    
    % Calculate residual and RMSE
    residual = (Y_val - Yb_hat_val);
    RMSE = rms(residual.^2);

end 