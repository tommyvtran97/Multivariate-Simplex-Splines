function [z_pred, XX_k1k1, dt] =  run_IEKF(U_k, Z_k, stdw, stdv, doIEKF)
    
    

    % Set simulation parameters
    n               = size(U_k, 2);
    dt              = 0.01;
    N               = 1000;
    epsilon         = 1e-10;
    maxIterations   = 100;

    %% Set Initial Values for States and Statistics
    Ex_0    = [Z_k(3,1); 0.5; 0.5; 0.5];  % initial estimate of optimal value of x_k_1k_1
    n       = length(stdw);         % number of states
    nm      = size(Z_k,1);          % number of measurements
    m       = size(U_k,2);          % number of inputs

    B       = eye(n);               % input matrix
    G       = eye(length(stdw));    % noise input matrix

    %% Initial Estimate for Covariance Matrix
    stdx_0  = [0.1, 0.1, 0.1, 0.1];     % Convergence KF depends on this estimate!
    P_0     = diag(stdx_0.^2);      % Create diagonal covariance matrix

    %% System Noise and Measurement Noise statistics
    Q = diag(stdw.^2);
    R = diag(stdv.^2);
    %w_k = stdw * randn(n, N) + Ew;  %TODO: Is this required?
    %v_k = stdv * randn(n, N) + Ev;  %TODO: Is this required?

    %% Initialize Arrays to Store Results
    XX_k1k1 = zeros(n, N);
    PP_k1k1 = zeros(n, N);
    STDx_cor = zeros(n, N);
    z_pred = zeros(nm, N);
    IEKFitcount = zeros(N, 1);

    x_k_1k_1 = Ex_0;            % x(0|0)=E{x_0}
    P_k_1k_1 = P_0;             % P(0|0)=P(0)

    % Extended Kalman Filter (EKF)
    ti = 0; 
    tf = dt;

    for k = 1:N
        % Prediction x(k+1|k) 
        [t, x_kk_1] = rk4(@kf_calc_f, x_k_1k_1,U_k(:,k), [ti tf]); 

        % Prediction Output z(k+1|k)
        z_kk_1 = kf_calc_h(0, x_kk_1, U_k(:,k)); %x_kk_1.^3; 
        z_pred(:,k) = z_kk_1;

        % Calc Phi(k+1,k) and Gamma(k+1, k)
        Fx = kf_calc_Fx(0, x_kk_1, U_k(:,k)); % perturbation of f(x,u,t)

        % The continuous to discrete time transformation of Df(x,u,t) and G
        [dummy, Psi] = c2d(Fx, B, dt);   
        [Phi, Gamma] = c2d(Fx, G, dt);   

        % Prediction covariance matrix P(k+1|k)
        P_kk_1 = Phi*P_k_1k_1*Phi' + Gamma*Q*Gamma'; 
        P_pred = diag(P_kk_1);
        stdx_pred = sqrt(diag(P_kk_1));

        % Run the Iterated Extended Kalman Filter
        if (doIEKF)

            % Iterative part
            eta2 = x_kk_1;
            err = 2*epsilon;

            itts = 0;
            while (err > epsilon)
                if (itts >= maxIterations)
                    fprintf("Terminating IEKF: exceeded max iterations (%d)\n", maxIterations);
                    break
                end
                itts = itts + 1;
                eta1 = eta2;

                Hx = kf_calc_Hx(0, eta1, U_k(:,k));

                % Check observability of state
                if (k == 1 && itts ==1)
                    rankHF = kf_calcObsRank(Hx,Fx);
                    if (rankHF < n)
                        warning('The current state is not observable: rank of Observability Matrix is %d, should be %d', rankHF, n);
                    end
                end

                % The innovation matrix
                Ve = (Hx*P_kk_1 * Hx' + R); 

                % Calculate the Kalman gain matrix
                K = P_kk_1 * Hx' / Ve;

                % New observation state
                z_p = kf_calc_h(0, eta1, U_k(:,k));

                eta2 = x_kk_1 + K * (Z_k(:,k) - z_p - Hx*(x_kk_1 - eta1));
                err  = norm((eta2 - eta1), inf) / norm(eta1, inf);

            end 

            IEFKitcount(k)  = itts;
            x_k_1k_1        = eta2;

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
        P_k_1k_1 = (eye(n) - K*Hx) * P_kk_1 * (eye(n) - K*Hx)' + K*R*K'; %[HUH?] 
        P_cor = diag(P_k_1k_1);
        stdx_cor = sqrt(diag(P_k_1k_1));

        % Next iteration
        ti = tf; 
        tf = tf + dt;

        % Store results
        XX_k1k1(:,k) = x_k_1k_1;
        %PP_k1k1(:,k) = P_k_1k_1;           % [Do I need this?]
        STDx_cor(:,k) = stdx_cor;
    end
    
end 