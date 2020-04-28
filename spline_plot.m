% Implement Simplex Spline Algorithm
clc;
clear;

simplex_order = 2;

num_triangles_x = 1;
num_triangles_y = 1;

% Define the grid boundaries
grid_begin_x    = -0.2;
grid_begin_y    = -0.2;
grid_end_x      = 0.8;
grid_end_y      = 0.2;

% Create Step Size
step_size_x     = (grid_end_x - grid_begin_x)/num_triangles_x;
step_size_y     = (grid_end_y - grid_begin_y)/num_triangles_y;

% Create Grid 
[x, y]          = meshgrid(grid_begin_x: step_size_x : grid_end_x,...
                grid_begin_y : step_size_y : grid_end_y);

Tri = delaunayTriangulation(x(:), y(:));
T = sort(Tri.ConnectivityList, 2)
multi_index     = sorted_bcoefficient(simplex_order);

plotID = 6001;
figure(plotID); 
hold on;
set(plotID, 'Position', [0 0 600 500], 'defaultaxesfontsize', 16, 'defaulttextfontsize', 14, 'color', [0.941, 0.941, 0.941], 'PaperPositionMode', 'auto');
set(gca, 'Color', [0.941, 0.941, 0.941]);
set(gca, 'XTick', [], 'YTIck', [],'XTickLabel',[],'YTickLabel',[]);
trimesh(T, x, y, [], 'EdgeColor', 'b', 'LineWidth', 2);
%title(sprintf('B-net (%d B-coefficients) for degree %d basis',size(multi_index,1), simplex_order), 'fontsize', 16)

vertices = Tri.Points;
for i = 1:size(vertices, 1)
    vertex_label = (['v_', num2str(i-0)]);
    text(vertices(i,1)+0.010, vertices(i,2), vertex_label, 'Color', 'red', 'FontSize', 15);
end

B_cart = [];
for i = 1:size(T, 1)
    BaryC = multi_index / simplex_order;
    simplex_coords  = vertices(T(i, :), :);
    B_cart          = vertcat([B_cart; bsplinen_bary2cart(simplex_coords, BaryC)]);
end

j = 0;
pos_x = -0.010
pos_y = -0.013
for i = 1:size(B_cart, 1)
    j = j + 1;
    plot(B_cart(i,1), B_cart(i,2), '.g', 'Markersize', 20)
    B_label = (['c_{' + string(multi_index(j,1)) + ',' + ...
        string(multi_index(j,2)) + ',' + string(multi_index(j,3)) + '}']);
    text(B_cart(i,1)-pos_x, B_cart(i,2)-pos_y, B_label, 'Color', 'black', 'FontSize', 12);
    
    if j == size(multi_index, 1)
        j = 0;
        pos_x = +0.070
        pos_y = +0.013
    end
end

