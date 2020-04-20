function IEKF_plot(alpha_m, beta_m, Vtot, z_pred, U_k, XX_k1k1, save)
    
    % Set simulation parameters
    n               = size(U_k, 2);
    dt              = 0.01;
    N               = 10000;
    epsilon         = 1e-10;
    doIEKF          = 1;
    maxIterations   = 100;
    
    state_u = XX_k1k1(1,:);
    state_v = XX_k1k1(2,:);
    state_w = XX_k1k1(3,:);
    state_C = XX_k1k1(4,:);
    
    kf_alpha = z_pred(1,:);
    kf_beta  = z_pred(2,:);
    kf_Vtot   = z_pred(3,:);
    
    % Construct the true angle attack by correcting for the bias
    alpha_true = z_pred(1,:) ./ (1 + state_C);
    
    % Define the time vector
    t = 0:dt:dt*size(U_k,2)-dt;

    % Show IEKF Results
    
    plotID = 1000;
    figure(plotID);
    subplot(311)
    plot(t(1:N), alpha_m(1:N), 'r-', 'LineWidth', 1.5)
    hold on
    plot(t(1:N), kf_alpha(1:N), 'b-', 'LineWidth', 1.5)
    hold on
    plot(t(1:N), alpha_true(1:N), 'g-', 'LineWidth', 1.5)
    legend('Measurement Data', 'Kalman Filter', 'True Data', 'location', 'northwest')
    xlim([0 50])
    set(gca,'FontSize',18);
    xlabel('time','interpreter','latex','FontSize', 25); 
    ylabel('$\alpha$ [rad]','interpreter','latex','FontSize', 25); 
    set(gcf, 'Position',  [0, 0, 1500, 1500])
    grid on

    subplot(312)
    plot(t(1:N), beta_m(1:N), 'r-', 'LineWidth', 1.5)
    hold on
    plot(t(1:N), kf_beta(1:N), 'b-', 'LineWidth', 1.5)
    xlim([0 50])
    set(gca,'FontSize',18);
    xlabel('time [s]','interpreter','latex','FontSize', 25); 
    ylabel('$\beta$ [rad]','interpreter','latex','FontSize', 25); 
    set(gcf, 'Position',  [0, 0, 1500, 1500])
    grid on
    
    subplot(313)
    plot(t(1:N), Vtot(1:N), 'r-', 'LineWidth', 1.5)
    hold on
    plot(t(1:N), kf_Vtot(1:N), 'b-', 'LineWidth', 1.5)
    xlim([0 50])
    set(gca,'FontSize',18);
    xlabel('time [s]','interpreter','latex','FontSize', 25); 
    label1 = ylabel('$V$ [m/s]','interpreter','latex','FontSize', 25); 
    set(gcf, 'Position',  [0, 0, 1500, 1500])
    label1.Position(1) = label1.Position(1) - 0.12;
    grid on
    if (save)
    saveas(gcf,[pwd,'\Plots\measurement_kf'],'epsc');
    end
    
    plotID = 1001;
    figure(plotID);
    subplot(411)
    plot(t(1:N), state_u(1:N), 'b-')
    xlim([0 50])
    set(gca,'FontSize',18);
    xlabel('time [s]','interpreter','latex','FontSize', 25); 
    ylabel('$u$ [m/s]','interpreter','latex','FontSize', 25); 
    set(gcf, 'Position',  [0, 0, 1500, 1500])
    grid on 

    subplot(412)
    plot(t(1:N), state_v(1:N), 'b-')
    xlim([0 50])
    set(gca,'FontSize',18);
    xlabel('time [s]','interpreter','latex','FontSize', 25); 
    label1 = ylabel('$v$ [m/s]','interpreter','latex','FontSize', 25); 
    set(gcf, 'Position',  [0, 0, 1500, 1500])
    label1.Position(1) = label1.Position(1) - 0.25;
    grid on
    
    subplot(413)
    plot(t(1:N), state_w(1:N), 'b-')
    xlim([0 50])
    set(gca,'FontSize',18);
    xlabel('time [s]','interpreter','latex','FontSize', 25); 
    label2 = ylabel('$w$ [m/s]','interpreter','latex','FontSize', 25); 
    set(gcf, 'Position',  [0, 0, 1500, 1500])
    label2.Position(1) = label2.Position(1) - 0.23;
    grid on

    subplot(414)
    plot(t(1:N), state_C(1:N), 'b-')
    xlim([0 50])
    set(gca,'FontSize',18);
    xlabel('time [s]','interpreter','latex','FontSize', 25); 
    label3 = ylabel('$C_{\alpha_{up}}$ [-]','interpreter','latex','FontSize', 25); 
    set(gcf, 'Position',  [0, 0, 1500, 1500])
    label3.Position(1) = label3.Position(1) - 0.21;
    grid on
    if (save)
    saveas(gcf,[pwd,'\Plots\state_kf'],'epsc');
    end
end
