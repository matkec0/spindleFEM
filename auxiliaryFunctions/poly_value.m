function [y] = poly_value(x, pol)
% returns a value of a polynomial pol (given by the array
% of its coefficients) at the point x

m = size(pol, 1);
n = size(pol, 2);
y = zeros(m, 1);

for i = 1 : n
    y = x * y + pol(:, n + 1 - i);
end

end