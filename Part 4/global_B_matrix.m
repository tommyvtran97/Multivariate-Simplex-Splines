% GLOBAL_B_MATRIX creates a regression matrix

function [global_B_id, global_B_val, global_idx_val, Y_hat_spline,...
    c_spline, VAR] = global_B_matrix(order, X_id, Y_id, X_val, Tri, T, H)
    
    % Initialize parameters
    X_id    = X_id';
    X_val   = X_val';
    Y_id    = Y_id';
    
    % Create B regression matrix
    [B_sorted] = sorted_bcoefficient(order);
    
    % Perform data membership search
    [IMapB_id, BaryB_id]    = tsearchn(Tri.Points, T, X_id);
    [IMapB_val, BaryB_val]  = tsearchn(Tri.Points, T, X_val);
    
    num_simplex = size(T, 1);
    
    % Initialize empty list to store results
    global_B_id     = [];
    global_B_val    = [];
    global_idx_id   = [];
    global_idx_val  = [];
    
    % Iterate over all simplices
    for k=1:1:num_simplex
        
        BaryB_id_new    = BaryB_id(find(IMapB_id == k), :);
        BaryB_val_new   = BaryB_val(find(IMapB_val == k), :);
        
        global_idx_id   = vertcat(global_idx_id, find(IMapB_id == k));
        global_idx_val  = vertcat(global_idx_val, find(IMapB_val == k));

        global_Bx_id    = (zeros(size(BaryB_id_new, 1), size(B_sorted, 1)));
        global_Bx_val   = (zeros(size(BaryB_val_new, 1), size(B_sorted, 1)));

        for i=1:1:size(global_Bx_id, 1)
            for j=1:1:size(global_Bx_id,2)
                % Define the multinomial
                multi_nom = factorial(order) / (factorial(B_sorted(j,1)) ...
                    * factorial(B_sorted(j,2)) * factorial(B_sorted(j,3)));

                global_Bx_id(i,j) = multi_nom * BaryB_id_new(i,1)^(B_sorted(j,1)) * ...
                    BaryB_id_new(i,2)^(B_sorted(j,2)) * BaryB_id_new(i,3)^(B_sorted(j,3));
            end
        end

        for i=1:1:size(global_Bx_val, 1)
            for j=1:1:size(global_Bx_val,2)
                % Define the multinomial
                multi_nom = factorial(order) / (factorial(B_sorted(j,1)) ...
                    * factorial(B_sorted(j,2)) * factorial(B_sorted(j,3)));

                global_Bx_val(i,j) = multi_nom * BaryB_val_new(i,1)^(B_sorted(j,1)) * ...
                    BaryB_val_new(i,2)^(B_sorted(j,2)) * BaryB_val_new(i,3)^(B_sorted(j,3));
            end
        end
        
        % Create diagonal block matrix
        global_B_id     = blkdiag(global_B_id, global_Bx_id);
        global_B_val    = blkdiag(global_B_val, global_Bx_val);
      
    end

    % Equality constraint OLS estimator
    matrix_A = pinv([global_B_id' * global_B_id, H'; ...
        H, zeros(size(H, 1), size(H,1))]);
    C1 = matrix_A(1:size(global_B_id,2), 1:size(global_B_id,2));
    c_spline = C1 * global_B_id' * Y_id(global_idx_id);
   
    % Calculate the variance
    COV = C1;   VAR = diag(C1);
   
    % Calculate estimated values from the simplex spline model
    Y_hat_spline = global_B_val * c_spline;

    % Check continuity condition
    epsilon = 10^(-4);
    if (max(H*c_spline)) < epsilon
        %fprintf('Continuity Satified!\n')
    else
        fprintf('Continuity Failed!\n')
    end
    
end