function [ret_value] = fixing(st, l)
% auxiliary function for graph of rods problem to get unique solution if there is no Dirichlet boundary condition

ret_value = poly_integral( polynomials(st, l), l );

end