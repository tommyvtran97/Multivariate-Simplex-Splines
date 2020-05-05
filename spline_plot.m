function [] = spline_plot(order, X_id, Y_id, X_val, Y_val,...
    Y_hat_spline, global_idx_val, plot_spline, save)
    
    if plot_spline
    
        % Initialize Parameters
        X_id    = X_id';
        X_val   = X_val';
        Y_id    = Y_id';
        Y_val   = Y_val';

        % Create triangulation
        Tri_spline = delaunayn(X_val(global_idx_val, [1 2]));
        Tri_val = delaunayn(X_val(:, [1 2]));

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

        % Set fancy options for plotting 
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
            figpath = 'Plots/';
            fpath = sprintf('3D_spline');
            savefname = strcat(figpath, fpath);
            print(plotID, '-dpng', '-r300', savefname);
        end
        
    end