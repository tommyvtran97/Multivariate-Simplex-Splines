function [] = validation_spline(order,spline_continuity, num_triangles_x, num_triangles_y, X_id, Y_id, X_val, Y_val, max_spline_order, max_continuity, max_simplices_xy, Y_hat_spline,...
    global_B_val, global_idx_val, c_spline, VAR, plot_spline,...
    plot_validation, save)

    if plot_validation
        % Fixed continuity and varying polynomial degree model performance
        [global_B_id, global_B_val, global_idx_val Y_hat_spline, c_spline, VAR, RMSE_x, RMSE_y] = simplex_continuity2(order, spline_continuity,...
            X_id, Y_id, X_val, Y_val, plot_spline);

        % Fixed polynomial degree and varying continuity model performance
        [global_B_id, global_B_val, global_idx_val Y_hat_spline, c_spline, VAR, RMSE_x_cont, RMSE_y_cont] = simplex_continuity3(order, max_continuity,...
            X_id, Y_id, X_val, Y_val, max_spline_order, plot_spline);

        %Varying number of simplices
        [global_B_id, global_B_val, global_idx_val Y_hat_spline, c_spline, VAR, RMSE_x_simp, RMSE_y_simp] = simplex_continuity4(order, spline_continuity,...
            max_simplices_xy, X_id, Y_id, X_val, Y_val, max_spline_order, plot_spline);
        
        [global_B_id, global_B_val, global_idx_val Y_hat_spline, c_spline, VAR, T, x, y, vertices, RMSE] = simplex_continuity1(order, spline_continuity,...
            num_triangles_x, num_triangles_y, X_id, Y_id, X_val, Y_val, plot_spline, save);
    end

    if plot_spline 
        
        % Initalize Parameters
        Y_val = Y_val';
        residual = Y_val(global_idx_val) - Y_hat_spline;

        conf = 1.96/sqrt(length(residual));
        [acx, lags] = xcorr(residual-mean(residual), 'coeff');

        coef_idx = 1:1:size(VAR, 1);
        coef_simplex_idx = 1:1:size(global_B_val, 2);
        
        if plot_validation
            plotID = 6003;
            figure(plotID);
            set(plotID, 'Position', [0 0 1500 500], 'defaultaxesfontsize', 16, 'defaulttextfontsize', 14, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
            subplot(131)
            plot(RMSE_x, RMSE_y, 'b^-', 'MarkerSize', 11, 'MarkerFaceColor',[0 0 1])
            ylabel('Root mean square error (RMSE) [-]','interpreter','latex');
            xlabel('Polynomial degree [-]','interpreter','latex');
            legend('Continuity 1', 'location', 'northwest','interpreter','latex');
            grid on;

            subplot(132)
            plot(RMSE_x_cont, RMSE_y_cont, 'b^-', 'MarkerSize', 11, 'MarkerFaceColor',[0 0 1])
            ylabel('Root mean square error (RMSE) [-]','interpreter','latex');
            xlabel('Continuity [-]','interpreter','latex');
            legend('Polynomial Degree 4', 'location', 'northwest','interpreter','latex');
            grid on;
            subplot(133)
            plot(RMSE_x_simp, RMSE_y_simp, 'b^-', 'MarkerSize', 11, 'MarkerFaceColor',[0 0 1])
            ylabel('Root mean square error (RMSE) [-]','interpreter','latex');
            xlabel('Number of simplices [-]','interpreter','latex');
            legend('Polynomial Degree 4 and Continuity 1', 'location', 'northwest','interpreter','latex');
            grid on;
            if (save)
                figpath = 'Plots/';
                fpath = sprintf('spline_rms_polycont');
                savefname = strcat(figpath, fpath);
                print(plotID, '-depsc', '-r300', savefname);
            end
        end

        plotID = 6004;
        figure(plotID);
        set(plotID, 'Position', [0 0 1500 500], 'defaultaxesfontsize', 16, 'defaulttextfontsize', 14, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
        subplot(121)
        nbins = 50;
        histogram(residual, nbins);
        ylabel('Number of residuals','interpreter','latex');
        xlabel('Residual $C_m$','interpreter','latex');
        legend(string(order) + 'th order spline', 'location', 'northwest');
        grid on;

        subplot(122)
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
            fpath = sprintf('Cm_auto_spline');
            savefname = strcat(figpath, fpath);
            print(plotID, '-dpng', '-r300', savefname);
        end

        plotID = 6005;
        figure(plotID);
        subplot(121)
        set(plotID, 'Position', [0 0 1500 500], 'defaultaxesfontsize', 16, 'defaulttextfontsize', 14, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
        bar(coef_simplex_idx, c_spline, 'b');
        xlabel('Index [-]', 'Interpreter', 'Latex');
        ylabel('Paramater coefficient [-]', 'Interpreter', 'Latex');
        legend(string(order)+ 'th order polynomial', 'location', 'northwest');
        grid on;

        subplot(122)
        bar(coef_idx, VAR, 'b');
        xlabel('Index [-]', 'Interpreter', 'Latex');
        ylabel('Coefficient variance', 'Interpreter', 'Latex');
        grid on;
        if (save)
            figpath = 'Plots/';
            fpath = sprintf('Spline_ParVar');
            savefname = strcat(figpath, fpath);
            print(plotID, '-depsc', '-r300', savefname);
        end
        
    end
    
end
    