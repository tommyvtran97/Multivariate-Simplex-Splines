function [Y_hat_id, Y_hat_val, theta_hat, A_matrix_val] = OLS_function(X_id, Y_id, X_val, expo)
    
    % Transpose matrices for correct for multiplication
    Y_id = Y_id';
    
    A_matrix_id  = zeros(size(X_id, 2), size(expo, 2));
    A_matrix_val = zeros(size(X_val, 2), size(expo, 2));
    
    for i=1:1:size(X_id,2)
        for j=1:1:size(expo,2)
            if (expo(1,j) ~= 0) & (expo(2,j) ~= 0)
                A_matrix_id(i,j) = X_id(1,i)^(expo(1,j))*X_id(2,i)^(expo(2,j));
            end
            if (expo(1,j) == 0)
                A_matrix_id(i,j) = X_id(2,i)^(expo(2,j));
            end
            if (expo(2,j) == 0)
                A_matrix_id(i,j) = X_id(1,i)^(expo(1,j));
            end
            if (expo(1,j) == 0) & (expo(2,j) == 0)
                A_matrix_id(i,j) = 1;
            end 
        end
    end 

    for i=1:1:size(X_val,2)
        for j=1:1:size(expo,2)
            if (expo(1,j) ~= 0) & (expo(2,j) ~= 0)
                A_matrix_val(i,j) = X_val(1,i)^(expo(1,j))*X_val(2,i)^(expo(2,j));
            end
            if (expo(1,j) == 0)
                A_matrix_val(i,j) = X_val(2,i)^(expo(2,j));
            end
            if (expo(2,j) == 0)
                A_matrix_val(i,j) = X_val(1,i)^(expo(1,j));
            end
            if (expo(1,j) == 0) & (expo(2,j) == 0)
                A_matrix_val(i,j) = 1;
            end 
        end
    end 
    
    % Using pinv to avoid matrix instability
    theta_hat = pinv(A_matrix_id) * Y_id;
    
    Y_hat_id = (A_matrix_id * theta_hat)';
    Y_hat_val = (A_matrix_val * theta_hat)';
    
end 