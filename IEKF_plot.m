function IEKF_plot(Z_k, Cm, z_pred, z_pred_correct, U_k, XX_k1k1,...
            IEKFitcount, plot_kalman, save)
    
    if (plot_kalman)
        
        % Initialize Parameters
        dt       = 0.01;
        N        = size(Z_k, 2);

        state_u  = XX_k1k1(1,:);
        state_v  = XX_k1k1(2,:);
        state_w  = XX_k1k1(3,:);
        state_C  = XX_k1k1(4,:);
        
        alpha_m  = Z_k(1,:); 
        beta_m   = Z_k(2,:);  
        Vtot     = Z_k(3,:);    
        
        kf_alpha = z_pred(1,:);
        kf_beta  = z_pred(2,:);
        kf_Vtot  = z_pred(3,:);

        % Construct the true angle attack by correcting for the bias
        alpha_true = z_pred_correct(1,:);

        % Define time vector
        t   = 0:dt:dt*size(U_k,2)-dt;
        idx = 1:1:size(IEKFitcount, 1);

        % Show IEKF Results
        plotID = 1001;
        figure(plotID);
        set(plotID, 'Position', [0 0 1500 1000], 'defaultaxesfontsize', 16, 'defaulttextfontsize', 14, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
        subplot(311)
        hold on;
        plot(t(1:N), alpha_m(1:N), 'r-', 'LineWidth', 1.5)
        plot(t(1:N), kf_alpha(1:N), 'b-', 'LineWidth', 1.5)
        plot(t(1:N), alpha_true(1:N), 'g-', 'LineWidth', 1.5)
        xlabel('time [s]','interpreter','latex'); 
        ylabel('$\alpha$ [rad]','interpreter','latex'); 
        legend('Measurement Data', 'Kalman Filter', 'True Data', 'location', 'northwest')
        grid on;
        
        subplot(312)
        hold on;
        plot(t(1:N), beta_m(1:N), 'r-', 'LineWidth', 1.5)
        plot(t(1:N), kf_beta(1:N), 'b-', 'LineWidth', 1.5)
        xlabel('time [s]','interpreter','latex'); 
        ylabel('$\beta$ [rad]','interpreter','latex'); 
        grid on;

        subplot(313)
        hold on; 
        plot(t(1:N), Vtot(1:N), 'r-', 'LineWidth', 1.5)
        plot(t(1:N), kf_Vtot(1:N), 'b-', 'LineWidth', 1.5)
        xlabel('time [s]','interpreter','latex'); 
        label1 = ylabel('$V$ [m/s]','interpreter','latex'); 
        label1.Position(1) = label1.Position(1) - 0.12;
        grid on;
        if (save)
        saveas(gcf,[pwd,'\Plots\KF_meas'],'epsc');
        end

        plotID = 1002;
        figure(plotID);
        set(plotID, 'Position', [0 0 1500 1000], 'defaultaxesfontsize', 16, 'defaulttextfontsize', 14, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
        subplot(411) 
        plot(t(1:N), state_u(1:N), 'b-')
        xlabel('time [s]','interpreter','latex'); 
        ylabel('$u$ [m/s]','interpreter','latex'); 
        grid on;

        subplot(412)
        plot(t(1:N), state_v(1:N), 'b-')
        xlabel('time [s]','interpreter','latex'); 
        label1 = ylabel('$v$ [m/s]','interpreter','latex'); 
        label1.Position(1) = label1.Position(1) - 0.25;
        grid on;

        subplot(413)
        plot(t(1:N), state_w(1:N), 'b-')
        xlabel('time [s]','interpreter','latex'); 
        label2 = ylabel('$w$ [m/s]','interpreter','latex'); 
        label2.Position(1) = label2.Position(1) - 0.23;
        grid on;

        subplot(414)
        plot(t(1:N), state_C(1:N), 'b-')
        xlabel('time [s]','interpreter','latex'); 
        label3 = ylabel('$C_{\alpha_{up}}$ [-]','interpreter','latex'); 
        label3.Position(1) = label3.Position(1) - 0.21;
        grid on;
        if (save)
        saveas(gcf,[pwd,'\Plots\KF_state'],'epsc');
        end

        plotID = 1003;
        figure(plotID); 
        set(plotID, 'Position', [0 0 600 500], 'defaultaxesfontsize', 16, 'defaulttextfontsize', 14, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
        hold on;
        plot3(alpha_m, beta_m, Cm, '-r');
        plot3(kf_alpha, kf_beta, Cm, '-b'); 
        plot3(alpha_true, kf_beta, Cm, '-g'); 
        ylabel('$\beta_m$ [rad]','interpreter','latex');
        xlabel('$\alpha_m$ [rad]','interpreter','latex');
        zlabel('C_m [-]');
        title('F16 CM(\alpha_m, \beta_m) raw datapoints only');
        legend('Measurement Data', 'Kalman Filter', 'True Data', 'location', 'northwest')
        view(0, 90);
        grid on; 
        if (save)
        saveas(gcf,[pwd,'\Plots\KF_3D'],'epsc');
        end

        plotID = 1004;
        figure(plotID);
        set(plotID, 'Position', [0 0 600 500], 'defaultaxesfontsize', 16, 'defaulttextfontsize', 14, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
        plot(idx, IEKFitcount, 'b-')
        xlabel('Data Index [-]','interpreter','latex');
        ylabel('Number of Iteration [-]','interpreter','latex');
        legend('Iterated Extended Kalman Filter', 'location', 'northwest');
        grid on; 
        if (save)
        saveas(gcf,[pwd,'\Plots\IEKF_iterations'],'epsc');
        end
    end 
end
