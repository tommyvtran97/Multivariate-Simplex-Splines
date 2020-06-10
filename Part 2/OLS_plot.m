% OLS_PLOT shows the results of polynomial obtained from the ordinary
% least square method
 
function [] = OLS_plot(order, X_val, Y_val, Y_hat_val, RMSE, plot_OLS, save)

    fprintf('Polynomial order %d and RMS: %.2f%%\n', order, RMSE);
    
    if (plot_OLS)
        
        % Initialize parameters
        Z_k         = X_val';
        alpha_m     = Z_k(:, 1);
        beta_m      = Z_k(:, 2);
        Cm          = Y_hat_val';
        Cm_val      = Y_val';
        residual    = (Y_val - Y_hat_val)';

        TRIeval = delaunayn(Z_k(:, [1 2]));

        % Show the results
        plotID = 2001;
        figure(plotID);
        set(plotID, 'Position', [0 0 1500 500], 'defaultaxesfontsize', 16, 'defaulttextfontsize', 14, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
        subplot(121)
        hold on;
        trisurf(TRIeval, alpha_m, beta_m, Cm, 'EdgeColor', 'none'); 
        plot3(alpha_m, beta_m, Cm_val, '.k'); 
        ylabel('beta [rad]');
        xlabel('alpha [rad]');
        zlabel('C_m [-]');
        title('F16 CM(\alpha_m, \beta_m) simple polynomial model of ' + string(order) + 'th order');
        view(140, 36);
        grid on;

        % Fancy options for plotting 
        set(gcf,'Renderer','OpenGL');
        poslight = light('Position',[0.5 .5 15],'Style','local');
        hlight = camlight('headlight');
        material([.3 .8 .9 25]);
        minz = min(Cm);
        shading interp;
        lighting phong;
        drawnow();

        subplot(122)
        hold on;
        trisurf(TRIeval, alpha_m, beta_m, Cm_val, 'EdgeColor', 'none'); 
        plot3(alpha_m, beta_m, Cm_val, '.k'); 
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
        minz = min(Cm_val);
        shading interp;
        lighting phong;
        drawnow();
        if (save)
            figpath = 'Plots/';
            fpath = sprintf('3D_OLS');
            savefname = strcat(figpath, fpath);
            print(plotID, '-dpng', '-r300', savefname);
        end 
        
    end

end 