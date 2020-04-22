function [VAR, COV] = OLS_plot(X_val, Y_val, Y_hat_val, MSE_x, MSE_y, A_matrix_val, theta_hat, save)
   
    % Initialize parameters
    Z_k = X_val';
    alpha_m = Z_k(:, 1);
    beta_m  = Z_k(:, 2);
    Cm = Y_hat_val';
    Cm_val = Y_val';
    residual = (Y_hat_val - Y_val)';
    
    TRIeval = delaunayn(Z_k(:, [1 2]));

    % Initialize viewing angle
    az = 140;
    el = 36;

    % Create figures
    plotID = 2001;
    figure(plotID);
    subplot(121)
    set(plotID, 'Position', [0 100 1500 500], 'defaultaxesfontsize', 18, 'defaulttextfontsize', 18, 'color', [1 1 1], 'PaperPositionMode', 'auto');
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
    
    plotID = 2002;
    figure(plotID);
    set(plotID, 'Position', [0 100 800 500], 'defaultaxesfontsize', 18, 'defaulttextfontsize', 18, 'color', [1 1 1], 'PaperPositionMode', 'auto');
    plot(MSE_x, MSE_y, 'b^-', 'MarkerSize', 14, 'MarkerFaceColor',[1 .6 .6])
    ylabel('Mean Squared Error (MSE) [-]');
    xlabel('Polynomial Degree [-]');
    legend('Validation Dataset Accuracy', 'location', 'northwest');
    grid on
    if (save)
    saveas(gcf,[pwd,'\Plots\OLS_accuracy'],'epsc');
    end

    plotID = 2003;
    figure(plotID);
    nbins = 50;
    set(plotID, 'Position', [0 100 800 500], 'defaultaxesfontsize', 18, 'defaulttextfontsize', 18, 'color', [1 1 1], 'PaperPositionMode', 'auto');
    histogram(residual, nbins)
    ylabel('Number of Residuals','interpreter','latex');
    xlabel('Residual $C_m$','interpreter','latex');
    legend('10th order polynomial', 'location', 'northwest');
    grid on
    if (save)
    saveas(gcf,[pwd,'\Plots\Cm_normal'],'epsc');
    end
    
    conf = 1.96/sqrt(length(residual));
    [acx, lags] = xcorr(residual-mean(residual), 'coeff');
    
    plotID = 2004;
    figure(plotID);
    set(plotID, 'Position', [0 100 800 500], 'defaultaxesfontsize', 18, 'defaulttextfontsize', 18, 'color', [1 1 1], 'PaperPositionMode', 'auto');
    line([lags(1), lags(end)], [conf, conf], 'Color','red','LineStyle','--')
    hold on
    line([lags(1), lags(end)], [-conf, -conf], 'Color','red','LineStyle','--')
    hold on 
    plot(lags, acx)
    xlabel('Number of lags', 'Interpreter', 'Latex')
    ylabel('Auto-correlation', 'Interpreter', 'Latex')
    legend('95% Confidence Interval')
    grid on
    if (save)
    saveas(gcf,[pwd,'\Plots\Cm_auto'],'epsc');
    end
    
    % The OLS Estimator is not BLUE since E{epsilon*epsilon'} ~= sigma*I
    COV = inv(A_matrix_val'*A_matrix_val)*A_matrix_val' * residual*residual' * A_matrix_val * inv(A_matrix_val'*A_matrix_val);
    VAR = diag(COV);
    
    coef_idx = 1:1:size(VAR, 1);
    coef_OLS_idx = 1:1:size(theta_hat, 1);
    
    plotID = 2005;
    figure(plotID);
    subplot(121)
    set(plotID, 'Position', [0 100 1500 500], 'defaultaxesfontsize', 18, 'defaulttextfontsize', 18, 'color', [1 1 1], 'PaperPositionMode', 'auto');
    bar(coef_OLS_idx, theta_hat, 'b')
    xlabel('Index [-]', 'Interpreter', 'Latex')
    ylabel('Paramater Coefficient [-]', 'Interpreter', 'Latex')
    legend('10th order polynomial', 'location', 'northwest')
    xlim([0 70])
    grid on
    
    subplot(122)
    bar(coef_idx, VAR, 'b')
    xlabel('Index [-]', 'Interpreter', 'Latex')
    ylabel('Coefficient Variance', 'Interpreter', 'Latex')
    xlim([0 70])
    grid on
    if (save)
    saveas(gcf,[pwd,'\Plots\OLS_ParVar'],'epsc');
    end

end 