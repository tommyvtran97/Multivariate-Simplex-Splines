function [] = validation(X_id, X_val, Y_id, Y_val, Y_hat_val, polynomial_order, max_polynomial_order, save, plot_OLS)
    
    if (plot_OLS)
 
        RMSE_x = [];
        RMSE_y = [];

        for order=1:1:max_polynomial_order
            expo = create_polynomial(order);
            [Y_hat_id, Y_hat_val, theta_hat, Ax_val] = OLS_function(X_id, Y_id, X_val, expo);

            residual  = (Y_val - Y_hat_val).^2;
            RMSE = rms(residual);

            RMSE_x = [RMSE_x, order];
            RMSE_y = [RMSE_y, RMSE];
        end
        
        residual = (Y_val - Y_hat_val)';

        expo = create_polynomial(polynomial_order);
        [Y_hat_id, Y_hat_val, theta_hat, Ax_val] = OLS_function(X_id, Y_id, X_val, expo);

        residual = (Y_val - Y_hat_val)';

        conf = 1.96/sqrt(length(residual));
        [acx, lags] = xcorr(residual-mean(residual), 'coeff');

        % Statistical Analysis
        sigma = (residual' * residual) / (size(Ax_val, 1) - size(Ax_val, 2));
        COV = sigma * pinv(Ax_val' * Ax_val);
        VAR = diag(COV);

        coef_idx = 1:1:size(VAR, 1);
        coef_OLS_idx = 1:1:size(theta_hat, 1);

        plotID = 3001;
        figure(plotID);
        set(plotID, 'Position', [0 100 800 500], 'defaultaxesfontsize', 18, 'defaulttextfontsize', 18, 'color', [1 1 1], 'PaperPositionMode', 'auto');
        plot(RMSE_x, RMSE_y, 'b^-', 'MarkerSize', 14, 'MarkerFaceColor',[1 .6 .6])
        ylabel('Root Mean Squared Error (RMSE) [-]');
        xlabel('Polynomial Degree [-]');
        legend('Validation Dataset Accuracy', 'location', 'northwest');
        grid on
        if (save)
        saveas(gcf,[pwd,'\Plots\OLS_accuracy'],'epsc');
        end

        plotID = 3002;
        figure(plotID);
        nbins = 50;
        set(plotID, 'Position', [0 100 800 500], 'defaultaxesfontsize', 18, 'defaulttextfontsize', 18, 'color', [1 1 1], 'PaperPositionMode', 'auto');
        histogram(residual, nbins)
        ylabel('Number of Residuals','interpreter','latex');
        xlabel('Residual $C_m$','interpreter','latex');
        legend(string(polynomial_order) + 'th order polynomial', 'location', 'northwest');
        grid on
        if (save)
        saveas(gcf,[pwd,'\Plots\Cm_normal'],'epsc');
        end

        plotID = 3003;
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

        plotID = 3004;
        figure(plotID);
        subplot(121)
        set(plotID, 'Position', [0 100 1500 500], 'defaultaxesfontsize', 18, 'defaulttextfontsize', 18, 'color', [1 1 1], 'PaperPositionMode', 'auto');
        bar(coef_OLS_idx, theta_hat, 'b')
        xlabel('Index [-]', 'Interpreter', 'Latex')
        ylabel('Paramater Coefficient [-]', 'Interpreter', 'Latex')
        legend(string(polynomial_order)+ 'th order polynomial', 'location', 'northwest')
        grid on

        subplot(122)
        bar(coef_idx, VAR, 'b')
        xlabel('Index [-]', 'Interpreter', 'Latex')
        ylabel('Coefficient Variance', 'Interpreter', 'Latex')
        grid on
        if (save)
        saveas(gcf,[pwd,'\Plots\OLS_ParVar'],'epsc');
        end
    end
    
end