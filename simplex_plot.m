function [] = simplex_plot(TRI, PHI, X_id, Y_id, X_val, Y_val, Yb_hat_val, simplex_order, plot_simplex, save)
    
    if (plot_simplex)
        
        X_val = X_val';
        Y_val = Y_val';
        
        plotID = 4001;
        figure(plotID);
        set(plotID, 'Position', [0 0 600 500], 'defaultaxesfontsize', 16, 'defaulttextfontsize', 14, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
        plot3(X_id(1,:), X_id(2,:), Y_id, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0 0], 'Markersize', 5);
        xlabel('$\alpha$','interpreter','latex');
        ylabel('$\beta$','interpreter','latex');
        zlabel('$C_m$','interpreter','latex');
        title(sprintf('F16 Model %d datapoints', size(X_id, 2)));
        view(20, 12);
        grid on;
        if (save)
        saveas(gcf,[pwd,'\Plots\3D_data'],'epsc');
        end
        
        plotID = 4002;
        figure(plotID);
        set(plotID, 'Position', [0 0 600 500], 'defaultaxesfontsize', 16, 'defaulttextfontsize', 14, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
        hold on;
        trimesh(TRI, PHI(:,1), PHI(:,2), 'Color', 'b', 'LineWidth', 3')
        plot(X_id(1,:), X_id(2,:), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0 0], 'Markersize', 5);
        xlabel('$\alpha$','interpreter','latex')
        ylabel('$\beta$','interpreter','latex')
        title('Simplex and Datapoints');
        grid on;
        if (save)
        saveas(gcf,[pwd,'\Plots\Simplex_data'],'epsc');
        end
        
        tri             = delaunayTriangulation(PHI);
        triangles       = sort(tri.ConnectivityList,2);
        multi_index     = sorted_bcoefficient(simplex_order);

        plotID = 4003;
        figure(plotID); 
        hold on;
        set(plotID, 'Position', [0 0 600 500], 'defaultaxesfontsize', 16, 'defaulttextfontsize', 14, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
        set(gca, 'Color', [0.941, 0.941, 0.941]);
        set(gca, 'XTick', [], 'YTIck', [],'XTickLabel',[],'YTickLabel',[]);
        trimesh(TRI, PHI(:,1), PHI(:,2), [], 'EdgeColor', 'b', 'LineWidth', 2);
        title(sprintf('B-net (%d B-coefficients) for degree %d basis',size(multi_index,1), simplex_order), 'fontsize', 16)

        vertices = PHI;
        for i = 1:size(vertices, 1)
            vertex_label = (['v_', num2str(i-0)]);
            text(vertices(i,1)-0.050, vertices(i,2)-0.013, vertex_label, 'Color', 'red', 'FontSize', 20);
        end 

        B_cart = [];
        for i = 1:size(vertices, 1)
            BaryC = multi_index / simplex_order;
            simplex_coord = vertices;
            B_cart = horzcat(bsplinen_bary2cart(simplex_coord, BaryC));
        end

        for i = 1:size(B_cart, 1)
            plot(B_cart(i,1), B_cart(i,2), '.g', 'Markersize', 20)
            B_label = (['c_{' + string(multi_index(i,1)) + ',' + ...
                string(multi_index(i,2)) + ',' + string(multi_index(i,3)) + '}']);
            text(B_cart(i,1), B_cart(i,2), B_label, 'Color', 'black', 'FontSize', 12);
        end
        if (save)
        saveas(gcf,[pwd,'\Plots\Bnet_simplex'],'epsc');
        end

        % Show 3D plots for the simplex polynomial model and validation model 
        
        % Initialize viewing angle
        az = 140;
        el = 36;

        % Create triangulation
        Tri_val     = delaunayn(X_val(:, [1 2]));
        
        plotID = 4004;
        figure(plotID);
        set(plotID, 'Position', [0 0 1500 500], 'defaultaxesfontsize', 16, 'defaulttextfontsize', 14, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
        subplot(121);
        hold on;
        trisurf(Tri_val, X_val(:,1), X_val(:,2), Yb_hat_val, 'EdgeColor', 'none'); 
        plot3(X_val(:,1), X_val(:,2), Y_val, '.k');
        view(az, el);
        ylabel('beta [rad]');
        xlabel('alpha [rad]');
        zlabel('C_m [-]');
        title('F16 CM(\alpha_m, \beta_m) Simplex Polynomial Order ' + string(simplex_order));
        view(az, el);
        grid on;

        % Set fancy options for plotting 
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
        view(az, el);
        grid on;

        % Set fancy options for plotting 
        set(gcf,'Renderer','OpenGL');
        poslight = light('Position',[0.5 .5 15],'Style','local');
        hlight = camlight('headlight');
        material([.3 .8 .9 25]);
        minz = min(Y_val);
        shading interp;
        lighting phong;
        drawnow();
        if (save)
        saveas(gcf,[pwd,'\Plots\3D_simplex'],'epsc');
        end
        
end