% RK4 calculates the one step-ahead prediction using Runge-Kutta method

function [t,w] = rk4(fn, xin, uin, t)
    
    % Initialize parameters
    g = 9.81; 
    a = t(1); 
    b = t(2);
    w = xin;
    N = 2;
    h = (b-a) / N;
    t = a;
    
    % Perform Runge-Kutta algorithm
    for j=1:N
        K1 = h * fn(t, w, uin);
        K2 = h * fn(t+h/2, w+K1/2, uin);
        K3 = h * fn(t+h/2, w+K2/2, uin);
        K4 = h * fn(t+h, w+K3, uin);

        w = w + (K1 + 2*K2 + 2*K3 + K4) / 6;
        t = a + j*h;
    end
    
end
