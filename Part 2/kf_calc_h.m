%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zpred = kf_calcHx(x) Calculates the output dynamics equation h(x,u,t) 
%   
% Adapted from C.C. de Visser Delft University of Technology 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function zpred = kf_calcHx(t, x, u)
    
    n        = size(x, 1);
    zpred    = zeros(3, 1);
    
    u = x(1); v = x(2); w = x(3); C = x(4);
    
    % Define the output dynamics 
    zpred(1) = atan(w/u) * (1 + C);         
    zpred(2) = atan(v/(sqrt(u^2 + w^2)));   
    zpred(3) = sqrt(u^2 + v^2 + w^2);      
end 
    