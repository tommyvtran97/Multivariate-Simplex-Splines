function [Y_hat_val, theta_hat, Ax_val] = OLS_function(X_id, Y_id, X_val, expo)
    
    % Transpose matrice
    Y_id = Y_id';
    
    % Initialize regression matrix
    Ax_id  = zeros(size(X_id, 2), size(expo, 2));
    Ax_val = zeros(size(X_val, 2), size(expo, 2));

    % Create the regression matrix
    for i=1:1:size(X_id,2)
        for j=1:1:size(expo,2)
            if (expo(1,j) ~= 0) & (expo(2,j) ~= 0)
                Ax_id(i,j) = X_id(1,i)^(expo(1,j))*X_id(2,i)^(expo(2,j));
            end
            if (expo(1,j) == 0)
                Ax_id(i,j) = X_id(2,i)^(expo(2,j));
            end
            if (expo(2,j) == 0)
                Ax_id(i,j) = X_id(1,i)^(expo(1,j));
            end
            if (expo(1,j) == 0) & (expo(2,j) == 0)
                Ax_id(i,j) = 1;
            end 
        end
    end 

    for i=1:1:size(X_val,2)
        for j=1:1:size(expo,2)
            if (expo(1,j) ~= 0) & (expo(2,j) ~= 0)
                Ax_val(i,j) = X_val(1,i)^(expo(1,j))*X_val(2,i)^(expo(2,j));
            end
            if (expo(1,j) == 0)
                Ax_val(i,j) = X_val(2,i)^(expo(2,j));
            end
            if (expo(2,j) == 0)
                Ax_val(i,j) = X_val(1,i)^(expo(1,j));
            end
            if (expo(1,j) == 0) & (expo(2,j) == 0)
                Ax_val(i,j) = 1;
            end 
        end
    end 
    
    % Calculate the parameter coefficients
    theta_hat = pinv(Ax_id' * Ax_id) * Ax_id' * Y_id;
    
    % Calculate the estimate values from OLS method
    Y_hat_val = (Ax_val * theta_hat)';   % Tranpose to correct size
    
end 