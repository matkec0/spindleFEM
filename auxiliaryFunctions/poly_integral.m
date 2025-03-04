function [int] = poly_integral(pol, bdry)
% returns the antiderivative polynomial of the polynomial "pol"
% with taking into account the border point "bdry"

[m, n] = size(pol);
int = zeros(m, 1);
aux = zeros(size(pol));

for i = 1 : n
    aux(:, i + 1) = pol(:, i) / i;  
end

for j = 1 : m
    int(j) = poly_value(bdry, aux(j, :));
end
 
end