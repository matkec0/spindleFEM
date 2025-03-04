function [rez] = poly_product(pol1, pol2)
% returns a polynomial which is equal to the product of
% polynomials pol1 and pol2

[m, n1] = size(pol1);
[~, n2] = size(pol2);
f = zeros(m, n1 + n2 - 1);

for i = 1 : n1
    for j = 1 : n2
        
        f(:, i + j - 1) = f(:, i + j - 1) + pol1(:, i) .* pol2(:, j);
    end
end

if m > 1
    rez = sum(f);
else
    rez = f;
end

end
