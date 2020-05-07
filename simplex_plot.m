% SIMPLEX_PLOT shows the results from the simplex polynomial model

function [] = simplex_plot(TRI, PHI, X_id, Y_id, X_val, Y_val, Yb_hat_val,...
    simplex_order, RMSE, plot_simplex, save)

    fprintf('Simplex order %d and RMS: %d\n', simplex_order, RMSE);
    
    if (plot_simplex)
        
        % Initialize parameters
        X_val           = X_val';
        Y_val           = Y_val';

        Tri_val         = delaunayn(X_val(:, [1 2]));
        multi_index     = sorted_bcoefficient(simplex_order);
        
        % Show the results
        plotID = 4001;
        figure(plotID);
        set(plotID, 'Position', [0 0 1500 500], 'defaultaxesfontsize', 16, 'defaulttextfontsize', 14, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
        subplot(121)
        hold on;
        set(gca, 'Color', [0.941, 0.941, 0.941]);
        set(gca, 'XTick', [], 'YTIck', [],'XTickLabel',[],'YTickLabel',[]);
        trimesh(TRI, PHI(:,1), PHI(:,2), [], 'EdgeColor', 'b', 'LineWidth', 2);
        title(sprintf('B-net (%d B-coefficients) for degree %d basis',size(multi_index,1), simplex_order), 'fontsize', 16)
        
        % Add labels for vertices
        vertices = PHI;
        for i = 1:size(vertices, 1)
            vertex_label = (['v_', num2str(i-0)]);
            text(vertices(i,1)+0.010, vertices(i,2), vertex_label, 'Color', 'red', 'FontSize', 15);
        end 

        B_cart = [];
        for i = 1:size(vertices, 1)
            BaryC = multi_index / simplex_order;
            simplex_coord = vertices;
            B_cart = horzcat(bsplinen_bary2cart(simplex_coord, BaryC));
        end
        
        % Add label for the B-coefficients
        for i = 1:size(B_cart, 1)
            plot(B_cart(i,1), B_cart(i,2), '.g', 'Markersize', 20)
            B_label = (['c_{' + string(multi_index(i,1)) + ',' + ...
                string(multi_index(i,2)) + ',' + string(multi_index(i,3)) + '}']);
            text(B_cart(i,1)-0.070, B_cart(i,2)-0.013, B_label, 'Color', 'black', 'FontSize', 12);
        end
        
        subplot(122)
        hold on;
        trimesh(TRI, PHI(:,1), PHI(:,2), 'Color', 'b', 'LineWidth', 3')
        plot(X_id(1,:), X_id(2,:), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0 0], 'Markersize', 5);
        xlabel('alpha [rad]','interpreter','latex')
        ylabel('beta [rad]','interpreter','latex')
        title('Simplex and F16 data');
        grid on; 
        if (save)
            figpath = 'Plots/';
            fpath = sprintf('Bnet_simplex_data');
            savefname = strcat(figpath, fpath);
            print(plotID, '-depsc', '-r300', savefname);
        end
        
        plotID = 4002;
        figure(plotID);
        set(plotID, 'Position', [0 0 1500 500], 'defaultaxesfontsize', 16, 'defaulttextfontsize', 14, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
        subplot(121);
        hold on;
        trisurf(Tri_val, X_val(:,1), X_val(:,2), Yb_hat_val, 'EdgeColor', 'none'); 
        plot3(X_val(:,1), X_val(:,2), Y_val, '.k');
        view(140, 36);
        ylabel('beta [rad]');
        xlabel('alpha [rad]');
        zlabel('C_m [-]');
        title('F16 CM(\alpha_m, \beta_m) simplex polynomial of ' + string(simplex_order) + 'th order');
        view(140, 36);
        grid on;

        % Fancy options for plotting 
        set(gcf,'Renderer','OpenGL');
        poslight = light('Position',[0.5 .5 15],'Style','local');
        hlight = camlight('headlight');
        material([.3 .8 .9 25]);
        minz = min(Yb_hat_val);
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
            fpath = sprintf('3D_simplex');
            savefname = strcat(figpath, fpath);
            print(plotID, '-dpng', '-r300', savefname);
        end
        
    end
    
end