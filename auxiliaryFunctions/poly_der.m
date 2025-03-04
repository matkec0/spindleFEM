function [der] = poly_der(f)
% returns a derivative of a polinomial f

[m, n] = size(f);
der = zeros(m, n - 1);

for i = 2 : n

    der(:, i - 1) = f(:, i) * (i - 1);
end

end