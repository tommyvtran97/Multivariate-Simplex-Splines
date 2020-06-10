% SORTED_BCOEFFICIENT creates a vector that consists of all possible
% permutations of a 3 digit number in lexographical order

function [B_sorted] = sorted_bcoefficient(simplex_order)
    
    % Initialize list
    B_sorted = [];
    
    % Create the algorithm
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
    
    % Transpose the matrix
    B_sorted = B_sorted';
    
end