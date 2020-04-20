function [Y_hat_id, Y_hat_val, Y_hat_full, theta_hat_id, theta_hat_val, theta_hat_full, A_matrix_id, A_matrix_val, A_matrix_full] = OLS_function(X_id, Y_id, X_val, Y_val, X_full, Y_full, expo)

    size(expo); % 2x10
    A_matrix_id     = zeros(size(X_id, 2), size(expo, 2));
    A_matrix_val    = zeros(size(X_val, 2), size(expo, 2));
    A_matrix_full   = zeros(size(X_full, 2), size(expo, 2));
    
    for i=1:1:size(X_id,2)
        for j=1:1:size(expo,2)
            if (expo(1,j) ~= 0) & (expo(2,j) ~= 0)
                A_matrix_id(i,j) = X_id(1,i)^(expo(1,j))*X_id(2,i)^(expo(2,j));
                A_matrix_val(i,j) = X_val(1,i)^(expo(1,j))*X_val(2,i)^(expo(2,j));
                A_matrix_full(i,j) = X_full(1,i)^(expo(1,j))*X_full(2,i)^(expo(2,j));
            end
            if (expo(1,j) == 0)
                A_matrix_id(i,j) = X_id(2,i)^(expo(2,j));
                A_matrix_val(i,j) = X_val(2,i)^(expo(2,j));
                A_matrix_full(i,j) = X_full(2,i)^(expo(2,j));
            end
            if (expo(2,j) == 0)
                A_matrix_id(i,j) = X_id(1,i)^(expo(1,j));
                A_matrix_val(i,j) = X_val(1,i)^(expo(1,j));
                A_matrix_full(i,j) = X_full(1,i)^(expo(1,j));
            end
            if (expo(1,j) == 0) & (expo(2,j) == 0)
                A_matrix_id(i,j) = 1;
                A_matrix_val(i,j) = 1;
                A_matrix_full(i,j) = 1;
            end 
        end
    end 
    
    theta_hat_id = inv(A_matrix_id' * A_matrix_id) * A_matrix_id' * Y_id;
    theta_hat_val = inv(A_matrix_val' * A_matrix_val) * A_matrix_val' * Y_val;
    theta_hat_full = inv(A_matrix_full' * A_matrix_full) * A_matrix_full' * Y_full;
    
    Y_hat_id = A_matrix_id * theta_hat_id;
    Y_hat_val = A_matrix_val * theta_hat_val;
    Y_hat_full = A_matrix_full * theta_hat_full;
end 