% CREATE_POLYNOMIAL creates a vector that contains all possible 2 digit
% permutations based on the order of the polynomial in a lexographical order

function [expo] = create_polynomial(poly_order)
    
    expo = [];
    for i=0:1:poly_order
        count = 0; 
        for j = 0:1:poly_order
            count = count + 1;
            if count == 1
                a = i;
                b = j;
                A = [a;b];
                expo = [expo, A];
            end 
            if count ~= 1
                a = a-1;
                b = j;
                A = [a;b];
                expo = [expo, A];
            end 
            if b == i
                break
            end 
        end 
    end 
    expo';
 
end 