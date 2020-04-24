function [] = simplex_plot(TRI, PHI, X_id, Y_id, X_val, Y_val, Yb_hat_val, simplex_order, plot_simplex, save)
    
    if (plot_simplex)
        X_val = X_val';
        Y_val = Y_val';
        
        plotID = 4001;
        figure(plotID);
        set(plotID, 'Position', [0 100 1000 1000], 'defaultaxesfontsize', 18, 'defaulttextfontsize', 18, 'color', [1 1 1], 'PaperPositionMode', 'auto');
        plot3(X_id(1,:), X_id(2,:), Y_id, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0 0], 'Markersize', 5);
        view(20, 12);
        grid on;
        xlabel('$\alpha$','interpreter','latex');
        ylabel('$\beta$','interpreter','latex');
        zlabel('$C_m$','interpreter','latex');
        title(sprintf('F16 Model %d datapoints', size(X_id, 2)));
        if (save)
        saveas(gcf,[pwd,'\Plots\3D_data'],'epsc');
        end

        plotID = 4002;
        figure(plotID);
        set(plotID, 'Position', [0 100 1000 1000], 'defaultaxesfontsize', 18, 'defaulttextfontsize', 18, 'color', [1 1 1], 'PaperPositionMode', 'auto');
        grid on;
        hold on;
        trimesh(TRI, PHI(:,1), PHI(:,2), 'Color', 'b', 'LineWidth', 3')
        plot(X_id(1,:), X_id(2,:), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0 0], 'Markersize', 5);
        xlabel('$\alpha$','interpreter','latex')
        ylabel('$\beta$','interpreter','latex')
        title('Simplex and Datapoints')
        if (save)
        saveas(gcf,[pwd,'\Plots\Simplex_data'],'epsc');
        end

        Tri = delaunayn(X_val(:, [1 2]));

        % Initialize viewing angle
        az = 140;
        el = 36;

        % Create triangulation
        Tri     = delaunayn(X_val(:, [1 2]));

        % Show 3D plots for the simplex polynomial model and validation model 
        plotID = 4003;
        figure(plotID);
        subplot(121);
        set(plotID, 'Position', [0 100 1500 500], 'defaultaxesfontsize', 18, 'defaulttextfontsize', 18, 'color', [1 1 1], 'PaperPositionMode', 'auto');
        trisurf(Tri, X_val(:,1), X_val(:,2), Yb_hat_val, 'EdgeColor', 'none'); 
        grid on;
        hold on;
        plot3(X_val(:,1), X_val(:,2), Y_val, '.k');
        view(az, el);
        ylabel('beta [rad]');
        xlabel('alpha [rad]');
        zlabel('C_m [-]');
        title('F16 CM(\alpha_m, \beta_m) Simplex Polynomial Order ' + string(simplex_order));

        % Set fancy options for plotting 
        set(gcf,'Renderer','OpenGL');
        hold on;
        poslight = light('Position',[0.5 .5 15],'Style','local');
        hlight = camlight('headlight');
        material([.3 .8 .9 25]);
        minz = min(Yb_hat_val);
        shading interp;
        lighting phong;
        drawnow();

        subplot(122)
        trisurf(Tri, X_val(:,1), X_val(:,2), Y_val, 'EdgeColor', 'none'); 
        grid on;
        hold on;
        plot3(X_val(:,1), X_val(:,2), Y_val, '.k');  
        view(az, el);
        ylabel('beta [rad]');
        xlabel('alpha [rad]');
        zlabel('C_m [-]');
        title('F16 CM(\alpha_m, \beta_m) Validation Model');

        % Set fancy options for plotting 
        set(gcf,'Renderer','OpenGL');
        hold on;
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