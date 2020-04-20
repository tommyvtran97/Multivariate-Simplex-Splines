function [] = OLS_plot(X_id, Y_id, Y_hat_id)
    
    size(X_id);
    size(Y_id);
   
    Z_k = X_id';
    
    alpha_m = Z_k(:, 1);
    beta_m  = Z_k(:, 2);
    Cm = Y_id';
    Cm1 = Y_hat_id;
    
    
    TRIeval = delaunayn(Z_k(:, [1 2]));

    %   viewing angles
    az = 140;
    el = 36;

    % create figures

    plotID = 201;
    figure(plotID);
    set(plotID, 'Position', [0 100 900 500], 'defaultaxesfontsize', 10, 'defaulttextfontsize', 10, 'color', [1 1 1], 'PaperPositionMode', 'auto');
    % plot data points
    plot3(alpha_m, beta_m, Cm, '.k'); % note that alpha_m = alpha, beta_m = beta, y = Cm
    view(0, 90); 
    ylabel('beta [rad]');
    xlabel('alpha [rad]');
    zlabel('C_m [-]');
    title('F16 CM(\alpha_m, \beta_m) raw datapoints only');
    
    plotID = 202;
    figure(plotID);
    set(plotID, 'Position', [0 100 900 500], 'defaultaxesfontsize', 10, 'defaulttextfontsize', 10, 'color', [1 1 1], 'PaperPositionMode', 'auto');
    % plot data points
    plot3(alpha_m, beta_m, Cm1, '.k'); % note that alpha_m = alpha, beta_m = beta, y = Cm
    view(0, 90); 
    ylabel('beta [rad]');
    xlabel('alpha [rad]');
    zlabel('C_m [-]');
    title('F16 CM(\alpha_m, \beta_m) raw datapoints only');
    
    plotID = 2001;
    figure(plotID);
    set(plotID, 'Position', [800 100 900 500], 'defaultaxesfontsize', 10, 'defaulttextfontsize', 10, 'color', [1 1 1], 'PaperPositionMode', 'auto');
    trisurf(TRIeval, alpha_m, beta_m, Cm, 'EdgeColor', 'none'); 
    grid on;
    hold on;
    % plot data points
    plot3(alpha_m, beta_m, Cm, '.k'); % note that alpha_m = alpha, beta_m = beta, y = Cm
    view(az, el);
    ylabel('beta [rad]');
    xlabel('alpha [rad]');
    zlabel('C_m [-]');
    title('F16 CM(\alpha_m, \beta_m) raw interpolation');
    % set fancy options for plotting 
    set(gcf,'Renderer','OpenGL');
    hold on;
    poslight = light('Position',[0.5 .5 15],'Style','local');
    hlight = camlight('headlight');
    material([.3 .8 .9 25]);
    minz = min(Cm);
    shading interp;
    lighting phong;
    drawnow();
    
    plotID = 2002;
    figure(plotID);
    set(plotID, 'Position', [800 100 900 500], 'defaultaxesfontsize', 10, 'defaulttextfontsize', 10, 'color', [1 1 1], 'PaperPositionMode', 'auto');
    trisurf(TRIeval, alpha_m, beta_m, Cm1, 'EdgeColor', 'none'); 
    grid on;
    hold on;
    % plot data points
    plot3(alpha_m, beta_m, Cm, '.k'); % note that alpha_m = alpha, beta_m = beta, y = Cm
    view(az, el);
    ylabel('beta [rad]');
    xlabel('alpha [rad]');
    zlabel('C_m [-]');
    title('F16 CM(\alpha_m, \beta_m) raw interpolation');
    % set fancy options for plotting 
    set(gcf,'Renderer','OpenGL');
    hold on;
    poslight = light('Position',[0.5 .5 15],'Style','local');
    hlight = camlight('headlight');
    material([.3 .8 .9 25]);
    minz = min(Cm);
    shading interp;
    lighting phong;
    drawnow();
    

end 