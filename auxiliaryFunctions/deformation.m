function [ret_value] = deformation(vertices, X)
% for given matrix "vertices" and solution vector "X" returns
% the deformation of the graph structure (sum of position and displacement);
% the solution is of the same form as matrix "vertices" (Nx3)

n = max(size(vertices,1));
ret_value = zeros(n, 3);

for i = 1 : n
    ret_value(i, :) = vertices(i, :) + X((i - 1) * 6 + 1 : (i - 1) * 6 + 3, 1)';
end

end