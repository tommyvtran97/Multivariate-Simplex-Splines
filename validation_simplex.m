function [] = validation_simplex(X_id, X_val, Y_id, Y_val, c_hat,...
    simplex_order, max_simplex_order, plot_result, save)
    
    if (plot_result)
 
        RMSE_x = [];
        RMSE_y = [];

        for order=1:1:max_simplex_order
            [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, RMSE] = simplex_polynomial(X_id, Y_id,...
                X_val, Y_val, order,...
                plot_result, save);
            RMSE_x = [RMSE_x, order];
            RMSE_y = [RMSE_y, RMSE];
        end

        [~, ~, ~, ~, Bx_val, ~, ~, ~, ~, residual, ~] = simplex_polynomial(X_id, Y_id,...
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

        plotID = 5001;
        figure(plotID);
        set(plotID, 'Position', [0 100 800 500], 'defaultaxesfontsize', 18, 'defaulttextfontsize', 18, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
        plot(RMSE_x, RMSE_y, 'b^-', 'MarkerSize', 14, 'MarkerFaceColor',[1 .6 .6])
        ylabel('Root Mean Squared Error (RMSE) [-]');
        xlabel('Polynomial Degree [-]');
        legend('Validation Dataset Accuracy', 'location', 'northwest');
        grid on
        if (save)
        saveas(gcf,[pwd,'\Plots\simplex_rms'],'epsc');
        end

        plotID = 5002;
        figure(plotID);
        nbins = 50;
        set(plotID, 'Position', [0 100 800 500], 'defaultaxesfontsize', 18, 'defaulttextfontsize', 18, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
        histogram(residual, nbins)
        ylabel('Number of Residuals','interpreter','latex');
        xlabel('Residual $C_m$','interpreter','latex');
        legend(string(simplex_order) + 'th order polynomial', 'location', 'northwest');
        grid on
        if (save)
        saveas(gcf,[pwd,'\Plots\Cm_normal'],'epsc');
        end

        plotID = 5003;
        figure(plotID);
        set(plotID, 'Position', [0 100 800 500], 'defaultaxesfontsize', 18, 'defaulttextfontsize', 18, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
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

        plotID = 5004;
        figure(plotID);
        subplot(121)
        set(plotID, 'Position', [0 100 1500 500], 'defaultaxesfontsize', 18, 'defaulttextfontsize', 18, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
        bar(coef_simplex_idx, c_hat, 'b')
        xlabel('Index [-]', 'Interpreter', 'Latex')
        ylabel('Paramater Coefficient [-]', 'Interpreter', 'Latex')
        legend(string(simplex_order)+ 'th order polynomial', 'location', 'northwest')
        grid on

        subplot(122)
        bar(coef_idx, VAR, 'b')
        xlabel('Index [-]', 'Interpreter', 'Latex')
        ylabel('Coefficient Variance', 'Interpreter', 'Latex')
        grid on
        if (save)
        saveas(gcf,[pwd,'\Plots\Simplex_ParVar'],'epsc');
        end
    end
    
end