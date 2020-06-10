% SPLINE_PLOT shows the results from the simplex spline model

function [] = spline_plot(order, continuity, X_id, Y_id, X_val, Y_val,...
    Y_hat_spline, global_idx_val, T, x, y, vertices, RMSE,...
    num_triangles_x, num_triangles_y, plot_spline, save)
    
    num_simplices = num_triangles_x*num_triangles_y*2;
    
    fprintf('Spline order %d, continuity %d, %d simplices and RMS: %.2f%%\n', order, continuity, num_simplices, RMSE);
    
    if plot_spline
        
        % Initialize parameters
        X_id    = X_id';
        X_val   = X_val';
        Y_id    = Y_id';
        Y_val   = Y_val';

        % Create triangulation
        Tri_spline      = delaunayn(X_val(global_idx_val, [1 2]));
        Tri_val         = delaunayn(X_val(:, [1 2]));
        
        % Calculate sorted B-coefficients
        multi_index     = sorted_bcoefficient(order);
        
        % Create index vector for all triangles
        vertex_index = [];
        for i=1:1:size(T,1)
            vertex_index = vertcat(vertex_index, multi_index);
        end
        
        % Show the results
        plotID = 6001;
        figure(plotID); 
        hold on;
        set(plotID, 'Position', [0 0 1500 1500], 'defaultaxesfontsize', 16, 'defaulttextfontsize', 14, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
        set(gca, 'Color', [0.941, 0.941, 0.941]);
        set(gca, 'XTick', [], 'YTIck', [],'XTickLabel',[],'YTickLabel',[]);
        trimesh(T, x, y, [], 'EdgeColor', 'b', 'LineWidth', 2);
        title(sprintf('B-net (%d B-coefficients per simplex) for degree %d basis',size(multi_index,1), order), 'fontsize', 16)
        
        
        % Sorting the triangulations
        T = sortrows(T,[1,-2]);
        
        % Create list with coordinates for the B-coefficient labels
        [label_position] = text_label(num_triangles_x, num_triangles_y, order);
         
       % Add vertex labels
        for i = 1:size(vertices, 1)
            vertex_label = (['v_{', num2str(i-0),'}']);
            text(vertices(i,1)+0.010, vertices(i,2), vertex_label, 'Color', 'red', 'FontSize', 15);
        end

        % Add B-coefficient labels
        B_cart = [];
        for i = 1:size(T, 1)
            BaryC = multi_index / order;
            simplex_coords  = vertices(T(i, :), :);
            B_cart          = vertcat([B_cart; bsplinen_bary2cart(simplex_coords, BaryC)]);
        end
        
        for i = 1:1:size(B_cart, 1)
            plot(B_cart(i,1), B_cart(i,2), '.g', 'Markersize', 20)
            B_label = (['c_{' + string(vertex_index(i,1)) + ',' + ...
                string(vertex_index(i,2)) + ',' + string(vertex_index(i,3)) + '}']);
            text(B_cart(i,1)+label_position(i,1), B_cart(i,2)+label_position(i,2), B_label, 'Color', 'black', 'FontSize', 12);
            
        end
            
        % Add triangle labels
        for k = 1:size(T,1)
            triangle_x = mean(vertices(T(k,:), 1));
            triangle_y = mean(vertices(T(k,:), 2));
            triangle_label = (['t_{', num2str(k-0), '}']);
            text(triangle_x, triangle_y, triangle_label, 'Color', 'black', 'FontSize', 15');
        end

        if save
            figpath = 'Plots/';
            fpath = sprintf('Spline_triangulation');
            savefname = strcat(figpath, fpath);
            print(plotID, '-depsc', '-r300', savefname);
        end

        plotID = 6002;
        figure(plotID);
        set(plotID, 'Position', [0 0 1500 500], 'defaultaxesfontsize', 16, 'defaulttextfontsize', 14, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
        subplot(121);
        hold on;
        trisurf(Tri_spline, X_val(global_idx_val,1), X_val(global_idx_val,2), Y_hat_spline, 'EdgeColor', 'none'); 
        plot3(X_val(global_idx_val,1), X_val(global_idx_val,2), Y_val(global_idx_val), '.k');
        ylabel('beta [rad]');
        xlabel('alpha [rad]');
        zlabel('C_m [-]');
        title('F16 CM(\alpha_m, \beta_m) simplex spline polynomial of ' + string(order) + 'th order');
        view(140, 36);
        grid on;

        % Fancy options for plotting 
        set(gcf,'Renderer','OpenGL');
        poslight = light('Position',[0.5 .5 15],'Style','local');
        hlight = camlight('headlight');
        material([.3 .8 .9 25]);
        minz = min(Y_hat_spline);
        shading interp;
        lighting phong;
        drawnow();

        subplot(122)
        hold on;
        trisurf(Tri_val, X_val(:,1), X_val(:,2), Y_val, 'EdgeColor', 'none'); 
        plot3(X_val(:,1), X_val(:,2), Y_val, '.k');  
        ylabel('beta [rad]');
        xlabel('alpha [rad]');
        zlabel('C_m [-]');
        title('F16 CM(\alpha_m, \beta_m) raw interpolation');
        view(140, 36);
        grid on;

        % Fancy options for plotting
        set(gcf,'Renderer','OpenGL');
        poslight = light('Position',[0.5 .5 15],'Style','local');
        hlight = camlight('headlight');
        material([.3 .8 .9 25]);
        minz = min(Y_val);
        shading interp;
        lighting phong;
        drawnow();
        if (save)
            figpath = 'Plots/';
            fpath = sprintf('3D_spline');
            savefname = strcat(figpath, fpath);
            print(plotID, '-dpng', '-r300', savefname);
        end
        
    end
    
end