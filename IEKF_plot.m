function IEKF_plot(alpha_m, beta_m, Vtot, z_pred, U_k, XX_k1k1)
    
    % Set simulation parameters
    n               = size(U_k, 2);
    dt              = 0.01;
    N               = 1000;
    epsilon         = 1e-10;
    doIEKF          = 1;
    maxIterations   = 100;

    kf_alpha = z_pred(1,:);
    kf_beta  = z_pred(2,:);
    kf_Vtot   = z_pred(3,:);
    
    state_u = XX_k1k1(1,:);
    state_v = XX_k1k1(2,:);
    state_w = XX_k1k1(3,:);
    state_C = XX_k1k1(4,:);

    t = 0:dt:dt*size(U_k,2)-dt;

    % Show IEKF Results
    figure(1)
    subplot(311)
    plot(t(1:N), alpha_m(1:N), 'g-')
    hold on
    plot(t(1:N), kf_alpha(1:N), 'b-')
    legend('Measurement', 'Prediction')
    xlim([0 10])

    subplot(312)
    plot(t(1:N), beta_m(1:N), 'g-')
    hold on
    plot(t(1:N), kf_beta(1:N), 'b-')
    xlim([0 10])

    subplot(313)
    plot(t(1:N), Vtot(1:N), 'g-')
    hold on
    plot(t(1:N), kf_Vtot(1:N), 'b-')
    xlim([0 10])
    
    figure(2)
    subplot(411)
    plot(t(1:N), state_u(1:N), 'g-')
    xlim([0 10])

    subplot(412)
    plot(t(1:N), state_v(1:N), 'g-')
    xlim([0 10])

    subplot(413)
    plot(t(1:N), state_w(1:N), 'g-')
    xlim([0 10])

    subplot(414)
    plot(t(1:N), state_C(1:N), 'g-')
    xlim([0 10])

end
