function [A] = a_forma(st, l)
% returns scalar elementary matrix for FEM for bilinear form
% a(u,v) = \int u' * v' dx,
% where l is the length of the domain, and
% st is the power of the polynomials ("k" in Pk elements)

A = zeros(st + 1, st + 1);
pol = polynomials(st, l);

for i = 1 : st + 1
    for j = 1 : st + 1       
        f1 = poly_der(pol(i, :));
        f2 = poly_der(pol(j, :));
        new = poly_product(f1, f2);
                
        A(i, j) = poly_integral(new, l);
    end
end

end