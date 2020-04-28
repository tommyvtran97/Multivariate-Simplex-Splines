function [global_] = global_B_matrix(order, X_id, X_val, Tri, T)
    
    % Initialize Parameters
    X_id = X_id';
    X_val = X_val';
    
    % Create B regression matrix
    [B_sorted] = sorted_bcoefficient(order);
    
    % Perform Member Search
    [IMapB_id, BaryB_id]  = tsearchn(Tri.Points, T, X_id);
    [IMapB_val, BaryB_val] = tsearchn(Tri.Points, T, X_val);
    
    num_simplex = size(T, 1);
   
    
    global_B = [];
    for k=1:1:num_simplex
        
        BaryB_id_new = BaryB_id(find(IMapB_id == k), :);
        BaryB_val_new = BaryB_val(find(IMapB_val == k), :);

        global_B_id      = (zeros(size(BaryB_id_new, 1), size(B_sorted, 1)));
        global_B_val     = (zeros(size(BaryB_val_new, 1), size(B_sorted, 1)));

        for i=1:1:size(global_B_id, 1)
            for j=1:1:size(global_B_id,2)
                % Define the multinomial
                multi_nom = factorial(size(B_sorted, 2)) / (factorial(B_sorted(j,1)) ...
                            * factorial(B_sorted(j,2)) * factorial(B_sorted(j,3)));

                global_B_id(i,j) = multi_nom * BaryB_id(i,1)^(B_sorted(j,1)) * ...
                    BaryB_id(i,2)^(B_sorted(j,2)) * BaryB_id(i,3)^(B_sorted(j,3));
            end
        end

        for i=1:1:size(global_B_val, 1)
            for j=1:1:size(global_B_val,2)
                % Define the multinomial
                multi_nom = factorial(size(B_sorted, 2)) / (factorial(B_sorted(j,1)) ...
                            * factorial(B_sorted(j,2)) * factorial(B_sorted(j,3)));

                global_B_val(i,j) = multi_nom * BaryB_val(i,1)^(B_sorted(j,1)) * ...
                    BaryB_val(i,2)^(B_sorted(j,2)) * BaryB_val(i,3)^(B_sorted(j,3));
            end
        end
        global_B = blkdiag(global_B, global_B_id);
        
    end

end