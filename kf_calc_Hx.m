%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H = kf_calcDHx(x) Calculates the Jacobian of the output dynamics equation f(x,u,t) 
%   
% Adapted from C.C. de Visser Delft University of Technology 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Hx = kf_calcDHx(t, x, u)

    u = x(1); v = x(2); w = x(3); C = x(4);

    % Calculate Jacobian matrix of output dynamics
    J11 = -w*(1+C)/(u^2 + w^2);     
    J12 = 0;                        
    J13 = u*(1+C)/(u^2 + w^2);      
    J14 = atan(w/u);                 
    
    J21 = -(u * v) / (sqrt(u^2 + w^2) * (u^2 + v^2 + w^2));
    J22 = sqrt(u^2 + w^2) / (u^2 + v^2 + w^2);
    J23 = -(w * v) / (sqrt(u^2 + w^2) * (u^2 + v^2 + w^2));
    J24 = 0;
    
    J31 = u / sqrt(u^2 + v^2 + w^2);
    J32 = v / sqrt(u^2 + v^2 + w^2);
    J33 = w / sqrt(u^2 + v^2 + w^2);
    J34 = 0;
    
    Hx  = [J11, J12, J13, J14;
            J21, J22, J23, J24;
            J31, J32, J33, J34];

end
