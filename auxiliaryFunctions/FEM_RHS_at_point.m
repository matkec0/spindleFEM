function [F] = FEM_RHS_at_point(fun, vr, points, N)
% returns a RHS for FEM system with forcing in one point of each bundle
% fun is a function handle
% vr are vertices
% points is array of indices of points to apply forcing
% N is size of FEM system (number of unknowns)

F = zeros(N,1);
for i = 1 : size(points,2)
    ind = 3*(points(i)-1)+1: 3*(points(i)-1)+3;
    F(ind) = fun(vr(points(i),:));
end
end

