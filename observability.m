function observability()

    syms('u', 'v', 'w','C_alpha_up', 'u_dot', 'v_dot', 'w_dot')

    f   = [u_dot; v_dot; w_dot; 0];
    h   = [atan(w/u) * (1 + C_alpha_up);         
            atan(v/(sqrt(u^2 + w^2)));   
            sqrt(u^2 + v^2 + w^2)];  
    x   = [u; v; w; C_alpha_up];
    x_0 = [150; 1; 1; 1];
    
    kf_calcNonlinObsRank(f, h, x, x_0);
    
end 
