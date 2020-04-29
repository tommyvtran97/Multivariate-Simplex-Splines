function [] = validation_spline(order, Y_hat_spline, Y_val,...
    global_B_val, global_idx_val, c_spline, VAR, plot_spline, save)
    
    if plot_spline 
        
        % Initalize Parameters
        Y_val = Y_val';
        residual = Y_val(global_idx_val) - Y_hat_spline;

        conf = 1.96/sqrt(length(residual));
        [acx, lags] = xcorr(residual-mean(residual), 'coeff');

        coef_idx = 1:1:size(VAR, 1);
        coef_simplex_idx = 1:1:size(global_B_val, 2);

    %     plotID = 5001;
    %     figure(plotID);
    %     set(plotID, 'Position', [0 0 600 500], 'defaultaxesfontsize', 16, 'defaulttextfontsize', 14, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
    %     plot(RMSE_x, RMSE_y, 'b^-', 'MarkerSize', 11, 'MarkerFaceColor',[0 0 1])
    %     ylabel('Root Mean Squared Error (RMSE) [-]','interpreter','latex');
    %     xlabel('Polynomial Degree [-]','interpreter','latex');
    %     legend('Validation Dataset Accuracy', 'location', 'northwest','interpreter','latex');
    %     grid on;
    %     if (save)
    %     saveas(gcf,[pwd,'\Plots\simplex_rms'],'epsc');
    %     end

        plotID = 6003;
        figure(plotID);
        set(plotID, 'Position', [0 0 600 500], 'defaultaxesfontsize', 16, 'defaulttextfontsize', 14, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
        nbins = 50;
        histogram(residual, nbins);
        ylabel('Number of Residuals','interpreter','latex');
        xlabel('Residual $C_m$','interpreter','latex');
        legend(string(order) + 'th order spline', 'location', 'northwest');
        grid on;
        if (save)
        saveas(gcf,[pwd,'\Plots\Cm_normal'],'epsc');
        end

        plotID = 6004;
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
        saveas(gcf,[pwd,'\Plots\Cm_auto'],'epsc');
        end

        plotID = 6005;
        figure(plotID);
        subplot(121)
        set(plotID, 'Position', [0 0 1500 500], 'defaultaxesfontsize', 16, 'defaulttextfontsize', 14, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
        bar(coef_simplex_idx, c_spline, 'b');
        xlabel('Index [-]', 'Interpreter', 'Latex');
        ylabel('Paramater Coefficient [-]', 'Interpreter', 'Latex');
        legend(string(order)+ 'th order polynomial', 'location', 'northwest');
        grid on;

        subplot(122)
        bar(coef_idx, VAR, 'b');
        xlabel('Index [-]', 'Interpreter', 'Latex');
        ylabel('Coefficient Variance', 'Interpreter', 'Latex');
        grid on;
        if (save)
        saveas(gcf,[pwd,'\Plots\Simplex_ParVar'],'epsc');
        end
        
    end
    
end
    