function [] = validation_simplex(X_id, X_val, Y_id, Y_val, c_hat,...
    simplex_order, max_simplex_order, plot_result, plot_validation, save)
    
    if (plot_result)
 
        if plot_validation
            RMSE_x = [];
            RMSE_y = [];

            for order=1:1:max_simplex_order
                [~, ~, ~, ~, ~, ~, ~, ~, RMSE] = simplex_polynomial(X_id, Y_id,...
                    X_val, Y_val, order,...
                    plot_result, save);
                RMSE_x = [RMSE_x, order];
                RMSE_y = [RMSE_y, RMSE];
            end
        end

        [~, ~, Bx_val, ~, ~, ~, ~, residual, ~] = simplex_polynomial(X_id, Y_id,...
            X_val, Y_val, simplex_order, plot_result, save);

        residual = residual';

        conf = 1.96/sqrt(length(residual));
        [acx, lags] = xcorr(residual-mean(residual), 'coeff');

        % Statistical Analysis
        sigma = (residual' * residual) / (size(Bx_val, 1) - size(Bx_val, 2));
        COV = sigma * pinv(Bx_val' * Bx_val);
        VAR = diag(COV);

        coef_idx = 1:1:size(VAR, 1);
        coef_simplex_idx = 1:1:size(Bx_val, 2);
        
        if plot_validation
            plotID = 5001;
            figure(plotID);
            set(plotID, 'Position', [0 0 1500 500], 'defaultaxesfontsize', 16, 'defaulttextfontsize', 14, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
            subplot(122)
            plot(RMSE_x, RMSE_y, 'b^-', 'MarkerSize', 11, 'MarkerFaceColor',[0 0 1])
            ylabel('Root mean squared error (RMSE) [-]','interpreter','latex');
            xlabel('Polynomial degree [-]','interpreter','latex');
            ylim([0.007 0.013]) 
            legend('Validation Dataset', 'location', 'northwest','interpreter','latex');
            grid on;

            subplot(121)
            nbins = 50;
            histogram(residual, nbins);
            ylabel('Number of residuals','interpreter','latex');
            xlabel('Residual $C_m$','interpreter','latex');
            legend(string(simplex_order) + 'th order polynomial', 'location', 'northwest');
            grid on;
            if (save)
                figpath = 'Plots/';
                fpath = sprintf('simplex_cm_rms');
                savefname = strcat(figpath, fpath);
                print(plotID, '-depsc', '-r300', savefname);
            end
        end
        
        plotID = 5002;
        figure(plotID);
        set(plotID, 'Position', [0 0 600 500], 'defaultaxesfontsize', 16, 'defaulttextfontsize', 14, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
        hold on;
        line([lags(1), lags(end)], [conf, conf], 'Color','red','LineStyle','--');
        line([lags(1), lags(end)], [-conf, -conf], 'Color','red','LineStyle','--');
        plot(lags, acx);
        xlabel('Number of lags', 'Interpreter', 'Latex');
        ylabel('Auto-correlation', 'Interpreter', 'Latex');
        legend('95% Confidence Interval');
        grid on;
        if (save)
            figpath = 'Plots/';
            fpath = sprintf('Cm_auto_simplex');
            savefname = strcat(figpath, fpath);
            print(plotID, '-depsc', '-r300', savefname);
        end

        plotID = 5003;
        figure(plotID);
        subplot(121)
        set(plotID, 'Position', [0 0 1500 500], 'defaultaxesfontsize', 16, 'defaulttextfontsize', 14, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
        bar(coef_simplex_idx, c_hat, 'b');
        xlabel('Index [-]', 'Interpreter', 'Latex');
        ylabel('Paramater coefficient [-]', 'Interpreter', 'Latex');
        legend(string(simplex_order)+ 'th order polynomial', 'location', 'northwest');
        grid on;

        subplot(122)
        bar(coef_idx, VAR, 'b');
        xlabel('Index [-]', 'Interpreter', 'Latex');
        ylabel('Coefficient variance', 'Interpreter', 'Latex');
        grid on;
        if (save)
            figpath = 'Plots/';
            fpath = sprintf('Simplex_ParVar');
            savefname = strcat(figpath, fpath);
            print(plotID, '-depsc', '-r300', savefname);
        end
    end
    
end