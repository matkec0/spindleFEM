function [X] = solutionsTransform (X1,vrN,edN,st1)
% permutes the values in the solution vector in a different way
% X1 is a solution of the FEM system
% vrN is number of vertices in graph
% edN is number of edges in graph
% st1 is a power (k in Pk elements) used for displacements

X = X1;
shift = 3*(vrN+ edN*(st1-1));
for ii = 1 : vrN + edN*(st1-1)
	X(6*(ii-1)+1:6*(ii-1)+3) = X1(3*(ii-1)+1:3*(ii-1)+3);
	X(6*(ii-1)+4:6*(ii-1)+6) = X1(shift+3*(ii-1)+1:shift+3*(ii-1)+3);
end
   
end