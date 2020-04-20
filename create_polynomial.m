function [exponential] = create_polynomial(N_states, poly_order)
    
    exponential = [];
    for i=0:1:poly_order
        for j = 0:1:poly_order
            if i+j > poly_order
                break
            end 
            A = [i;j];
            exponential = [exponential, A];
        end 
    end 
    exponential';
 
end 