function [H, Tri, T] = simplex_continuity(order, continuity, num_triangles_x, num_triangles_y)

    % Define the grid boundaries
    grid_begin_x    = -0.2;
    grid_begin_y    = -0.2;
    grid_end_x      = 0.8;
    grid_end_y      = 0.2;

    % Create Step Size
    step_size_x     = (grid_end_x - grid_begin_x)/num_triangles_x;
    step_size_y     = (grid_end_y - grid_begin_y)/num_triangles_y;

    % Create Grid 
    [x, y]          = meshgrid(grid_begin_x: step_size_x : grid_end_x,...
                    grid_begin_y : step_size_y : grid_end_y);

    Tri             = delaunayTriangulation(x(:), y(:));
    T               = sort(Tri.ConnectivityList, 2);
    multi_index     = sorted_bcoefficient(order);
    
    % Create index vector for all triangles
    vertex_index = [];
    for i=1:1:size(T,1)
        vertex_index = vertcat(vertex_index, multi_index);
    end

    plotID = 6001;
    figure(plotID); 
    hold on;
    set(plotID, 'Position', [0 0 600 500], 'defaultaxesfontsize', 16, 'defaulttextfontsize', 14, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
    set(gca, 'Color', [0.941, 0.941, 0.941]);
    set(gca, 'XTick', [], 'YTIck', [],'XTickLabel',[],'YTickLabel',[]);
    trimesh(T, x, y, [], 'EdgeColor', 'b', 'LineWidth', 2);
    
    % Add vertex labels
    vertices = Tri.Points;
   
    for i = 1:size(vertices, 1)
        vertex_label = (['v_', num2str(i-0)]);
        text(vertices(i,1)+0.010, vertices(i,2), vertex_label, 'Color', 'red', 'FontSize', 15);
    end
 
    % Add B-coefficient labels
    B_cart = [];
    for i = 1:size(T, 1)
        BaryC = multi_index / order;
        simplex_coords  = vertices(T(i, :), :);
        B_cart          = vertcat([B_cart; bsplinen_bary2cart(simplex_coords, BaryC)]);
    end
    
    j = 0;
    pos_x = +0.070;
    pos_y = +0.013;
    for i = 1:1:size(B_cart, 1)
        j = j + 1;
        plot(B_cart(i,1), B_cart(i,2), '.g', 'Markersize', 20)
        B_label = (['c_{' + string(vertex_index(j,1)) + ',' + ...
            string(vertex_index(j,2)) + ',' + string(vertex_index(j,3)) + '}']);
        text(B_cart(i,1)-pos_x, B_cart(i,2)-pos_y, B_label, 'Color', 'black', 'FontSize', 12);
    end
    
        
    %Add triangle labels
    for k = 1:size(T,1)
        triangle_x = mean(vertices(T(k,:), 1));
        triangle_y = mean(vertices(T(k,:), 2));
        triangle_label = (['t_', num2str(k-0)]);
        text(triangle_x, triangle_y, triangle_label, 'Color', 'black', 'FontSize', 15');
    end
%--------------------------------------------------------------------------------------------    
    
    % Find all the edges for continuity
    int_edges = setdiff(sort(edges(Tri),2), sort(freeBoundary(Tri),2), 'rows');
    
    % Find triangle with interior edge
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

   % Implement smoothness matrix
   H = [];
   
   for i=1:1:size(int_edges, 1)
       size_LH = size(multi_index, 1) * (triangle_edge_list(i, 1)-1);
       size_RH = size(multi_index, 1) * (triangle_edge_list(i, 2)-1);
       for order=0:1:continuity
           
            % Permutations left hand side
            multi_index_LH = multi_index(find(multi_index(:,index_OOE(i,1)) == order), :);
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
            gamma = sorted_bcoefficient(order);
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
           
            % Create smoothness matrix
            smooth_matrix = zeros(size(multi_index_LH, 1), size(multi_index, 1)*size(T,1));

            OOE_BaryV1 = OOE_vertex_bary_1(i,:);

            % Store result into smoothness matrix
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
          








