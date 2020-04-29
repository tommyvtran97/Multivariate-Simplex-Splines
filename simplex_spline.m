function [H, Tri, T] = simplex_spline(order, continuity, num_triangles_x, num_triangles_y)

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
    T               = sort(sort(Tri.ConnectivityList, 2), 1);
    multi_index     = sorted_bcoefficient(order);

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
        B_label = (['c_{' + string(multi_index(j,1)) + ',' + ...
            string(multi_index(j,2)) + ',' + string(multi_index(j,3)) + '}']);
        text(B_cart(i,1)-pos_x, B_cart(i,2)-pos_y, B_label, 'Color', 'black', 'FontSize', 12);

        if j == size(multi_index, 1)
            j = 0;
            pos_x = -0.010;
            pos_y = -0.013;
        end
    end
    
    % Implement Smoothness Matrix 
    %% TODO: high order continuity & multi simplices
    H = [];
    
    for order=0:1:continuity

        % Permutations left hand side
        multi_index_LH = multi_index(find(multi_index(:,1) == order), :);

        % Permutation right hand side
        gamma = sorted_bcoefficient(order);
        zero_vector = zeros(size(multi_index_LH, 1), 1);
        matrix_RH  = horzcat(multi_index_LH(:,2), multi_index_LH(:,3), zero_vector);

        multi_index_RH = [];
        for i=1:1:size(matrix_RH, 1)
            for j=1:1:size(gamma, 1)
                multi_index_RH = vertcat(multi_index_RH, matrix_RH(i,:) + gamma(j,:));
            end
        end 

        % Create smoothness matrix
        smooth_matrix = zeros(size(multi_index_LH, 1), size(multi_index, 1)*2);
        
        OOE_V1 = Tri.Points(T(1,1), :);
        OOE_BaryV1 = bsplinen_cart2bary(Tri.Points(T(2,:),:), OOE_V1);
        
        % Store result into smoothness matrix
        index = 1;
        k = 1;
        for i=1:1:size(multi_index_LH, 1)
            [~, idx_left] = ismember(multi_index_LH(i,:), multi_index, 'rows');
            smooth_matrix(i, idx_left) = -1;

            for j=index:1:size(multi_index_RH, 1)
                    [~, idx_right] = ismember(multi_index_RH(j,:), multi_index, 'rows');
                    if k <= size(OOE_BaryV1, 2)
                        smooth_matrix(i, idx_right+size(multi_index, 1)) = OOE_BaryV1(1, k);
                    end

                    if k == size(gamma, 1)
                        index = j + 1;
                        k = 1;
                        break
                    end
                    k = k + 1;
            end
        end
        % Concatenate smoothness matrix vertically
        H = vertcat(H, smooth_matrix);
    end
end








