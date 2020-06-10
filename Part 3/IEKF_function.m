% IEKF_FUNCTION performs the Iterated Extended Kalman Filter

function [z_pred, z_pred_corr, XX_k1k1, IEKFitcount, n] = ...
            IEKF_function(U_k, Z_k, stdw, stdv, IEKF)
        
    % Set simulation parameters
    dt              = 0.01;
    N               = size(U_k, 2);
    epsilon         = 1e-10;
    maxIterations   = 100;

    % Set initial values for states and statistics
    Ex_0    = [Z_k(3,1); 0.5; 0.5; 0.5];    % initial estimate of optimal value of x_k_1k_1
    n       = length(stdw);                 % number of states
    nm      = size(Z_k,1);                  % number of measurements
    m       = size(U_k,1);                  % number of inputs

    B       = eye(n);                       % input matrix
    G       = eye(length(stdw));            % noise input matrix

    % Initial estimate for covariance matrix
    stdx_0  = [sqrt(0.1), sqrt(0.1), sqrt(0.1), sqrt(0.1)];     % Convergence KF depends on this estimate!
    P_0     = diag(stdx_0.^2);                                  % Create diagonal covariance matrix
    
    %% System Noise and Measurement Noise statistics
    Q = diag(stdw.^2);
    R = diag(stdv.^2);

    %% Initialize Arrays to Store Results
    XX_k1k1     = zeros(n, N);
    PP_k1k1     = zeros(n, N);
    STDx_cor    = zeros(n, N);
    z_pred      = zeros(nm, N);
    IEKFitcount = zeros(N, 1);

    x_k_1k_1 = Ex_0;            % x(0|0)=E{x_0}
    P_k_1k_1 = P_0;             % P(0|0)=P(0)

    % Extended Kalman Filter (EKF)
    ti = 0; 
    tf = dt;

    for k = 1:N
        % Prediction x(k+1|k) 
        [t, x_kk_1] = rk4(@kf_calc_f, x_k_1k_1, U_k(:,k), [ti tf]);

        % Prediction uutput z(k+1|k)
        z_kk_1      = kf_calc_h(0, x_kk_1, U_k(:,k));
        z_pred(:,k) = z_kk_1;

        % Calculate Phi(k+1,k) and Gamma(k+1, k)
        Fx = kf_calc_Fx(0, x_kk_1, U_k(:,k)); % perturbation of f(x,u,t)

        % The continuous to discrete time transformation of Df(x,u,t) and G  
        [Phi, Gamma] = c2d(Fx, G, dt);

        % Prediction covariance matrix P(k+1|k)
        P_kk_1      = Phi * P_k_1k_1 * Phi' + Gamma * Q * Gamma';
        P_pred      = diag(P_kk_1);
        stdx_pred   = sqrt(diag(P_kk_1));

        % Run the Iterated Extended Kalman Filter
        if (IEKF)

            % Iterative part
            eta2 = x_kk_1;
            err  = 2 * epsilon;

            itts = 0;
            while (err > epsilon)
                if (itts >= maxIterations)
                    fprintf("Terminating IEKF: exceeded max iterations (%d)\n",...
                    maxIterations)
                    break
                end
                
                itts = itts + 1;
                eta1 = eta2;

                Hx = kf_calc_Hx(0, eta1, U_k(:,k));

                % The innovation matrix
                Ve = (Hx * P_kk_1 * Hx' + R);

                % Calculate the Kalman gain matrix
                K = P_kk_1 * Hx' / Ve;

                % New observation state
                z_p = kf_calc_h(0, eta1, U_k(:,k));

                eta2 = x_kk_1 + K * (Z_k(:,k) - z_p - Hx*(x_kk_1 - eta1));
                err  = norm((eta2 - eta1), inf) / norm(eta1, inf);

            end 

            IEKFitcount(k,1)  = itts;
            x_k_1k_1          = eta2;

        else  
            % Correction
            Hx = kf_calc_Hx(0, x_kk_1, U_k(:,k)); 

            % Pz(k+1|k) (covariance matrix of innovation)
            Ve = (Hx*P_kk_1 * Hx' + R); 

            % K(k+1) (gain)
            K = P_kk_1 * Hx' / Ve;

            % Calculate optimal state x(k+1|k+1) 
            x_k_1k_1 = x_kk_1 + K * (Z_k(:,k) - z_kk_1); 

        end

        % Calculate covariance matrix of state estimation error
        P_k_1k_1    = (eye(n) - K*Hx) * P_kk_1 * (eye(n) - K*Hx)' + K*R*K'; 
        P_cor       = diag(P_k_1k_1);
        stdx_cor    = sqrt(diag(P_k_1k_1));

        % Next iteration
        ti = tf; 
        tf = tf + dt;

        % Store results
        XX_k1k1(:,k)   = x_k_1k_1;          
        STDx_cor(:,k)  = stdx_cor;
    end
    
    % Correct for bias alpha
    z_pred_corr = z_pred;
    z_pred_corr(1,:) = z_pred_corr(1,:) ./ (1 + XX_k1k1(4, end));
    
end 