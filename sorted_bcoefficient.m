function [B_sorted] = sorted_coefficient(simplex_order)
   
    B_sorted = [];
    % Write algorithm
    for i=simplex_order:-1:0
        for k = 0:1:simplex_order
            for j = simplex_order:-1:0
                if (i+j+k) == simplex_order
                    A = [i;j;k];
                    B_sorted = [B_sorted, A];
                end
            end
        end
    end
    B_sorted = B_sorted';
end