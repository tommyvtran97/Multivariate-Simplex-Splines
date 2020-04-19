%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F = kf_calcDFx(x) Calculates the Jacobian of the system dynamics equation f(x,u,t) 
%   
%   Adapted from C.C. de Visser Delft University of Technology 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DFx = kf_calcDFx(t, x, u)

    % Calculate Jacobian matrix of system dynamics
    DFx = zeros(4);      % Derivative is zero hence 4x4 matrix with zeros        
    
end
    

