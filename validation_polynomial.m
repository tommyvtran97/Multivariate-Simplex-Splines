% VALIDATION_POLYNOMIAL shows the results of the validation of the
% polynomial model.

function [] = validation_polynomial(X_id, Y_id, X_val, Y_val,...
    polynomial_order, max_polynomial_order, plot_result, plot_validation, save)
    
    if (plot_result)
        
        if plot_validation
            RMSE_x = [];
            RMSE_y = [];

            for order=1:1:max_polynomial_order
                [Y_hat_val, theta_hat, Ax_val, ~] = OLS_function(order,X_id,...
                    Y_id, X_val, Y_val);

                residual  = (Y_val - Y_hat_val);
                RMSE = rms(residual);

                RMSE_x = [RMSE_x, order];
                RMSE_y = [RMSE_y, RMSE];
            end
            
        end
        
        % Calculate the residual
        [Y_hat_val, theta_hat, Ax_val, ~] = OLS_function(polynomial_order, X_id, Y_id, X_val,Y_val);
        residual = (Y_val - Y_hat_val)';
        
        % Calculate the autocorrelation of the residual
        conf = 1.96/sqrt(length(residual));
        [acx, lags] = xcorr(residual-mean(residual), 'coeff');

        % Statistical analysis
        M1 = pinv(Ax_val' * Ax_val) * Ax_val';
        M2 = Ax_val * pinv(Ax_val' * Ax_val);
        E_residual = cov(residual);
        COV = M1*E_residual*M2;
        VAR = diag(COV);

        coef_idx = 1:1:size(VAR, 1);
        coef_OLS_idx = 1:1:size(theta_hat, 1);
        
        if plot_validation
            
            plotID = 3001;
            figure(plotID);
            set(plotID, 'Position', [0 0 1500 500], 'defaultaxesfontsize', 16, 'defaulttextfontsize', 14, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
            subplot(121)
            nbins = 50;
            histogram(residual, nbins);
            xlabel('Residual $C_m$','interpreter','latex');
            ylabel('Number of residuals','interpreter','latex');
            legend(string(polynomial_order) + 'th order polynomial', 'location', 'northwest');
            grid on;

            subplot(122)
            plot(RMSE_x, RMSE_y, 'b^-', 'MarkerSize', 11, 'MarkerFaceColor',[0 0 1]);
            xlabel('Polynomial order [-]','interpreter','latex');
            ylabel('Root mean square error (RMSE) [-]','interpreter','latex');
            ylim([0.007 0.016]);
            legend('Validation Dataset', 'location', 'northwest');
            grid on;
            if (save)
                figpath = 'Plots/';
                fpath = sprintf('OLS_residual_rms');
                savefname = strcat(figpath, fpath);
                print(plotID, '-dpng', '-r300', savefname);
            end
            
        end
        
        plotID = 3002;
        figure(plotID);
        set(plotID, 'Position', [0 0 1500 500], 'defaultaxesfontsize', 16, 'defaulttextfontsize', 14, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
        subplot(121)
        hold on;
        line([lags(1), lags(end)], [conf, conf], 'Color','red','LineStyle','--');
        line([lags(1), lags(end)], [-conf, -conf], 'Color','red','LineStyle','--');
        plot(lags, acx);
        xlabel('Number of lags', 'Interpreter', 'Latex');
        ylabel('Auto-correlation', 'Interpreter', 'Latex');
        legend('95% Confidence Interval');
        grid on;
        
        subplot(122)
        plot(Y_val, residual, '.k')
        ylabel('Residual [-]', 'Interpreter', 'Latex');
        xlabel('Predicted value [-]', 'Interpreter', 'Latex');
        grid on;
        if (save)
            figpath = 'Plots/';
            fpath = sprintf('Cm_auto_simple_variance');
            savefname = strcat(figpath, fpath);
            print(plotID, '-depsc', '-r300', savefname);
        end

        plotID = 3003;
        figure(plotID);
        subplot(121)
        set(plotID, 'Position', [0 0 1500 500], 'defaultaxesfontsize', 16, 'defaulttextfontsize', 14, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
        bar(coef_OLS_idx, theta_hat, 'b');
        xlabel('Index [-]', 'Interpreter', 'Latex');
        ylabel('Paramater coefficient [-]', 'Interpreter', 'Latex');
        legend(string(polynomial_order)+ 'th order polynomial', 'location', 'northwest');
        grid on;

        subplot(122)
        bar(coef_idx, VAR, 'b');
        xlabel('Index [-]', 'Interpreter', 'Latex');
        ylabel('Coefficient variance', 'Interpreter', 'Latex');
        grid on;
        if (save)
            figpath = 'Plots/';
            fpath = sprintf('OLS_ParVar');
            savefname = strcat(figpath, fpath);
            print(plotID, '-depsc', '-r300', savefname);
        end
        
    end
    
end