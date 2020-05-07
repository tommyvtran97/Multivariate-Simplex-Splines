% SIMPLEX_CONTINUITY1 creates a smoothness matrix that is required for 
% the continuity conditions and also creates the global regression matrix.
% Furthermore, it iterates over different continuity orders to perform
% a quality analysis of the model.

function [global_B_id, global_B_val, global_idx_val Y_hat_spline, ...
    c_spline, VAR, RMSE_x_cont, RMSE_y_cont] = simplex_continuity(...
    spline_order, max_continuity, X_id, Y_id, X_val, Y_val,...
    max_spline_order, plot_spline)
    
    % Create empty list for RMS for increasing continuity order
    RMSE_x_cont = [];
    RMSE_y_cont = [];
    
    % Define the number of simplices in x and y direction
    num_triangles_x = 1;
    num_triangles_y = 1;
    
    for max_order_continuity=0:1:max_continuity
        
        % Define the grid boundaries
        grid_begin_x    = -0.2;
        grid_begin_y    = -0.2;
        grid_end_x      = 0.8;
        grid_end_y      = 0.2;

        % Create step size of the grid
        step_size_x     = (grid_end_x - grid_begin_x)/num_triangles_x;
        step_size_y     = (grid_end_y - grid_begin_y)/num_triangles_y;

        % Create the grid
        [x, y]          = meshgrid(grid_begin_x: step_size_x : grid_end_x,...
                        grid_begin_y : step_size_y : grid_end_y);

        Tri             = delaunayTriangulation(x(:), y(:));
        T               = sort(Tri.ConnectivityList, 2);
        multi_index     = sorted_bcoefficient(max_spline_order);
        vertices        = Tri.Points;

        % Create index vector for all triangles
        vertex_index = [];
        for i=1:1:size(T,1)
            vertex_index = vertcat(vertex_index, multi_index);
        end

        % Find all the edges for each simplex
        int_edges = setdiff(sort(edges(Tri),2), sort(freeBoundary(Tri),2), 'rows');

        % Find the matching simplices for each edge
        triangle_edge_list = [];
        triangle_list = [];
        for i=1:1:size(int_edges, 1)
            triangle_count = 0;
            for j=1:1:size(T, 1)
                vertex_count = 0;
                for k=1:1:size(T, 2)
                    if T(j,k) == int_edges(i,1) || T(j,k) == int_edges(i, 2)
                        vertex_count = vertex_count + 1;
                        if vertex_count == 2
                            triangle_count = triangle_count + 1;
                            triangle_list = [triangle_list, j];
                            if length(triangle_list) == 2
                                triangle_edge_list = vertcat(triangle_edge_list, triangle_list);
                                triangle_list = [];
                            end
                        end
                    end
                end
            end
        end

        % Calculate for each set of triangles the out of edge vertex
        OOE_vertex_bary_1 = [];
        OOE_vertex_bary_2 = [];
        for i=1:1:size(triangle_edge_list, 1)
           triangle_1 = vertices(T(triangle_edge_list(i, 1), :), :);
           triangle_2 = vertices(T(triangle_edge_list(i, 2), :), :);

           edge = vertices(int_edges(i,:), :);

           OOE_cart_1 = setdiff(triangle_1,  edge, 'rows');
           OOE_cart_2 = setdiff(triangle_2,  edge, 'rows');

           OOE_bary_1 = bsplinen_cart2bary(triangle_2, OOE_cart_1);
           OOE_bary_2 = bsplinen_cart2bary(triangle_1, OOE_cart_2);

           OOE_vertex_bary_1 = vertcat(OOE_vertex_bary_1, OOE_bary_1);
           OOE_vertex_bary_2 = vertcat(OOE_vertex_bary_2, OOE_bary_2);
        end 

        % Calculate correct out of edge index to use Boor's continuity equation
        index_OOE = [];
        index_list = [];
        for i=1:1:size(triangle_edge_list, 1)
           for j=1:1:size(triangle_edge_list, 2)
               for k =1:1:size(OOE_vertex_bary_1, 2)
                   if T(triangle_edge_list(i, j), k) ~= int_edges(i,1) && T(triangle_edge_list(i, j), k) ~= int_edges(i,2)
                       index_list = [index_list, k];
                   end
               end
           end
           index_OOE = vertcat(index_OOE, index_list);
           index_list = [];
        end

        % Implement the smoothness matrix
        H = [];

        for i=1:1:size(int_edges, 1)
           size_LH = size(multi_index, 1) * (triangle_edge_list(i, 1)-1);
           size_RH = size(multi_index, 1) * (triangle_edge_list(i, 2)-1);
           for order_continuity=0:1:max_order_continuity

                % Permutations left hand side
                multi_index_LH = multi_index(find(multi_index(:,index_OOE(i,1)) == order_continuity), :);
                if index_OOE(i,1) == 1
                    kappa_0 = multi_index_LH(:, 2);
                    kappa_1 = multi_index_LH(:, 3);
                end
                if index_OOE(i,1) == 2
                    kappa_0 = multi_index_LH(:, 1);
                    kappa_1 = multi_index_LH(:, 3);
                end
                if index_OOE(i,1) == 3
                    kappa_0 = multi_index_LH(:, 1);
                    kappa_1 = multi_index_LH(:, 2);
                end

                % Permutation right hand side
                gamma = sorted_bcoefficient(order_continuity);
                if index_OOE(i,2) == 1
                    column_1 = zeros(size(multi_index_LH, 1), 1);
                    column_2 = kappa_0;
                    column_3 = kappa_1;
                end
                if index_OOE(i,2) == 2
                    column_1 = kappa_0;
                    column_2 = zeros(size(multi_index_LH, 1), 1);
                    column_3 = kappa_1;
                end
                if index_OOE(i,2) == 3
                    column_1 = kappa_0;
                    column_2 = kappa_1;
                    column_3 = zeros(size(multi_index_LH, 1), 1);
                end


                matrix_RH  = horzcat(column_1, column_2, column_3);

                multi_index_RH = [];
                for k=1:1:size(matrix_RH, 1)
                    for j=1:1:size(gamma, 1)
                        multi_index_RH = vertcat(multi_index_RH, matrix_RH(k,:) + gamma(j,:));
                    end
                end

                % Create the smoothness matrix
                smooth_matrix = zeros(size(multi_index_LH, 1), size(multi_index, 1)*size(T,1));

                OOE_BaryV1 = OOE_vertex_bary_1(i,:);

                % Store result into the smoothness matrix
                index = 1;
                id = 1;
                for k=1:1:size(multi_index_LH, 1)
                    [~, idx_left] = ismember(multi_index_LH(k,:), multi_index, 'rows');
                    smooth_matrix(k, idx_left+size_LH) = -1;

                    for j=index:1:size(multi_index_RH, 1)
                            [~, idx_right] = ismember(multi_index_RH(j,:), multi_index, 'rows');
                            if id <= size(gamma,1)
                                smooth_matrix(k, idx_right+size_RH) =...
                                    OOE_BaryV1(1, 1)^gamma(id, 1)*OOE_BaryV1(1, 2)^gamma(id, 2)*...
                                    OOE_BaryV1(1, 3)^gamma(id, 3);
                            end

                            if id == size(gamma, 1)
                               index = j + 1;
                               id = 1;
                               break
                            else
                               id = id + 1;
                            end
                    end
                end
                % Concatenate smoothness matrix vertically
                H = vertcat(H, smooth_matrix);

           end
        end

        % Create global B regression matrix
        [global_B_id, global_B_val, global_idx_val, Y_hat_spline, c_spline, VAR]...
        = global_B_matrix(max_spline_order, X_id, Y_id, X_val, Tri, T, H);
        
        % Calculate the root mean square (RMS)
        residual = Y_val(global_idx_val)' - Y_hat_spline;
        RMSE = rms(residual);

        RMSE_x_cont = [RMSE_x_cont, order_continuity];
        RMSE_y_cont = [RMSE_y_cont, RMSE];
    end
end
