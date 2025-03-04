function [A] = b2_forma(st1, st2, l)
% returns scalar elementary matrix for FEM for bilinear form
% b1(u,v) = \int u * v dx,
% where l is the length of the domain, and
% st1 and st2 are powers of the polynomials for u and v
% respectively ("k" in Pk elements)


A = zeros(st1 + 1, st2 + 1);
pol1 = polynomials(st1, l);
pol2 = polynomials(st2, l);

for i = 1 : st1 + 1
    for j = 1 : st2 + 1
        f1 = pol1(i, :);
        f2 = pol2(j, :);
        new = poly_product(f1, f2);
                
        A(i, j) = poly_integral(new, l);
    end
end

end