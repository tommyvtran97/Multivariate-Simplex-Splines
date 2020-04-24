function [] = OLS_plot(X_val, Y_val, Y_hat_val, save, plot_OLS)
    
    if (plot_OLS)
        % Initialize parameters
        Z_k = X_val';
        alpha_m = Z_k(:, 1);
        beta_m  = Z_k(:, 2);
        Cm = Y_hat_val';
        Cm_val = Y_val';
        residual = (Y_val - Y_hat_val)';

        TRIeval = delaunayn(Z_k(:, [1 2]));

        % Initialize viewing angle
        az = 140;
        el = 36;

        % Create figures
        plotID = 2001;
        figure(plotID);
        subplot(121)
        set(plotID, 'Position', [0 100 1500 500], 'defaultaxesfontsize', 18, 'defaulttextfontsize', 18, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
        trisurf(TRIeval, alpha_m, beta_m, Cm, 'EdgeColor', 'none'); 
        grid on;
        hold on;
        plot3(alpha_m, beta_m, Cm_val, '.k'); 
        view(az, el);
        ylabel('beta [rad]');
        xlabel('alpha [rad]');
        zlabel('C_m [-]');
        title('F16 CM(\alpha_m, \beta_m) OLS Simple Polynomial Model');

        % Set fancy options for plotting 
        set(gcf,'Renderer','OpenGL');
        hold on;
        poslight = light('Position',[0.5 .5 15],'Style','local');
        hlight = camlight('headlight');
        material([.3 .8 .9 25]);
        minz = min(Cm);
        shading interp;
        lighting phong;
        drawnow();

        subplot(122)
        trisurf(TRIeval, alpha_m, beta_m, Cm_val, 'EdgeColor', 'none'); 
        grid on;
        hold on;
        plot3(alpha_m, beta_m, Cm_val, '.k'); 
        view(az, el);
        ylabel('beta [rad]');
        xlabel('alpha [rad]');
        zlabel('C_m [-]');
        title('F16 CM(\alpha_m, \beta_m) raw interpolation');

        % Set fancy options for plotting 
        set(gcf,'Renderer','OpenGL');
        hold on;
        poslight = light('Position',[0.5 .5 15],'Style','local');
        hlight = camlight('headlight');
        material([.3 .8 .9 25]);
        minz = min(Cm_val);
        shading interp;
        lighting phong;
        drawnow();
        if (save)
        saveas(gcf,[pwd,'\Plots\OLS_3D'],'epsc');
        end 
    end
    

end 