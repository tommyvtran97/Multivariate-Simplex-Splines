function IEKF_plot(Z_k, z_pred, z_pred_correct, XX_k1k1,...
            IEKFitcount, plot_kalman, save)
    
    if (plot_kalman)
        
        % Initialize parameters
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
        t   = 0:dt:dt*size(Z_k,2)-dt;
        idx = 1:1:size(IEKFitcount, 1);

        % Show IEKF results
        
        plotID = 1001;
        figure(plotID);
        set(plotID, 'Position', [0 0 1500 600], 'defaultaxesfontsize', 18, 'defaulttextfontsize', 14, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
        subplot(311)
        hold on;
        plot(t(1:N), alpha_m(1:N), 'r-', 'LineWidth', 1.5)
        plot(t(1:N), kf_alpha(1:N), 'b-', 'LineWidth', 1.5)
        plot(t(1:N), alpha_true(1:N), 'g-', 'LineWidth', 1.5)
        xlabel('Time [s]','interpreter','latex'); 
        label1 = ylabel('$\alpha$ [rad]','interpreter','latex','FontSize', 20); 
        label_ref_1 = label1.Position(1) - 0.2;
        legend('Measurement Data', 'Iterated Extended Kalman Filter', 'True Data', 'location', 'northwest')
        grid on;
        
        subplot(312)
        hold on;
        plot(t(1:N), beta_m(1:N), 'r-', 'LineWidth', 1.5)
        plot(t(1:N), kf_beta(1:N), 'b-', 'LineWidth', 1.5)
        xlabel('Time [s]','interpreter','latex','FontSize', 20); 
        label1 = ylabel('$\beta$ [rad]','interpreter','latex','FontSize', 20);
        label1.Position(1) = label_ref_1 - 0.2;
        grid on;

        subplot(313)
        hold on; 
        plot(t(1:N), Vtot(1:N), 'r-', 'LineWidth', 1.5)
        plot(t(1:N), kf_Vtot(1:N), 'b-', 'LineWidth', 1.5)
        xlabel('Time [s]','interpreter','latex','FontSize', 20); 
        label1 = ylabel('$V$ [m/s]','interpreter','latex','FontSize', 20); 
        label1.Position(1) = label_ref_1 - 0.2;
        grid on;
        if (save)
            figpath = 'Plots/';
            fpath = sprintf('KF_meas');
            savefname = strcat(figpath, fpath);
            print(plotID, '-depsc', '-r300', savefname);
        end

        plotID = 1002;
        figure(plotID);
        set(plotID, 'Position', [0 0 1500 600], 'defaultaxesfontsize', 18, 'defaulttextfontsize', 16, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
        subplot(411) 
        plot(t(1:N), state_u(1:N), 'b-')
        xlabel('Time [s]','interpreter','latex','FontSize', 20); 
        label1 = ylabel('$u$ [m/s]','interpreter','latex','FontSize', 20); 
        label_ref_1 = label1.Position(1);
        grid on;

        subplot(412)
        plot(t(1:N), state_v(1:N), 'b-')
        xlabel('Time [s]','interpreter','latex','FontSize', 20); 
        label1 = ylabel('$v$ [m/s]','interpreter','latex','FontSize', 20); 
        label1.Position(1) = label_ref_1;
        grid on;

        subplot(413)
        plot(t(1:N), state_w(1:N), 'b-')
        xlabel('Time [s]','interpreter','latex','FontSize', 20); 
        label1 = ylabel('$w$ [m/s]','interpreter','latex','FontSize', 20); 
        label1.Position(1) = label_ref_1;
        grid on;

        subplot(414)
        plot(t(1:N), state_C(1:N), 'b-')
        xlabel('Time [s]','interpreter','latex','FontSize', 20); 
        label1 = ylabel('$C_{\alpha_{up}}$ [-]','interpreter','latex','FontSize', 20); 
        label1.Position(1) = label_ref_1;
        grid on;
        if (save)
            figpath = 'Plots/';
            fpath = sprintf('KF_state');
            savefname = strcat(figpath, fpath);
            print(plotID, '-depsc', '-r300', savefname);
        end

        plotID = 1003;
        figure(plotID);
        set(plotID, 'Position', [0 0 600 500], 'defaultaxesfontsize', 18, 'defaulttextfontsize', 14, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
        plot(idx, IEKFitcount, 'b-')
        xlabel('Data index [-]','interpreter','latex');
        ylabel('Number of iteration [-]','interpreter','latex');
        legend('Iterated Extended Kalman Filter', 'location', 'northwest');
        grid on; 
        if (save)
            figpath = 'Plots/';
            fpath = sprintf('IEKF_iterations');
            savefname = strcat(figpath, fpath);
            print(plotID, '-depsc', '-r300', savefname);
        end
    end 
end
