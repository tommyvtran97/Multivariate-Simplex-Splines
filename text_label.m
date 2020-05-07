% TEXT_LABEL creates a list with coordinates for the text labels of the
% B-coefficient which is used for the triagulation plots


function [pos_label] = text_label(num_triangles_x, num_triangles_y, spline_order)
    
    % Initialize position
    x = 0.01;
    y = 0.01;
    
    % Initialize scaling paramaters
    sx = 3;
    sxx = 3;
    sy = 3;
    syy = 2;
    
    % Create position matrix to the size of spline order
    M1_pos1 = [];   M1_pos2 = [];
    M2_pos1 = [];   M2_pos2 = [];

    for i=1:1:spline_order
        M1 = zeros(i,2);
        M1(i,2) = y;
        M1_pos1 = vertcat(M1_pos1, M1);
        
        if i ~= 1
            M1 = zeros(i,2);
            M1(1,2) = y;
            M1(end, 2) = -y;
            M1_pos2 = vertcat(M1_pos2, M1);
        end
    end
    
    for i=1:1:spline_order
        M2 = zeros(1,2);
        M2(1,2) = -y;;
        M2_pos1 = vertcat(M2_pos1, M2);
        
        if i ~= 1
            M2 = zeros(1,2);
            M2(1,1) = -sx*x;
            M2(1,2) = 0;
            M2_pos2 = vertcat(M2_pos2, M2);
        end
    end
    
    label_position_1 = [M1_pos1; M2_pos1; -sxx*x y;];
    label_position_2 = [sx*x -y; M1_pos2; -sx*x syy*y; M2_pos2; -sx*x -y;];
    label_position_comb = vertcat(label_position_1, label_position_2);
    
    pos_label = [];
    for i=1:1:num_triangles_y
        pos_label = vertcat(pos_label, label_position_comb);
    end
    
    pos_temp_label = pos_label;
    
    for k=1:1:num_triangles_x-1
        pos_label = vertcat(pos_label, pos_temp_label);
    end

end