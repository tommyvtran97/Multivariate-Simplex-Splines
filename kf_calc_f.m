%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xdot = kf_calc_Fx(x, u ,t) Calculates the system dynamics equation
% f(x,u,t) 
%   
% Adapted from C.C. de Visser Delft University of Technology 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xdot = kf_calcFx(t, x, u)

    n = size(x, 1);
    xdot = zeros(n, 1);
    
    % system dynamics go here!
    xdot = [u; 0];
end
