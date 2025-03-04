function [pol] = polynomials(st, l)
% for a power "st" of the FEM elements (the k in Pk) and
% the length l of the interval, returns a matrix
% (st+1) x (st+1) filled with the coefficients of all 
% basis polynomials

if st == 0
    pol = 1;
    return;
end

pol = ones(st + 1, st + 1);

for i = 1 : st + 1
    aux = 1;
    fac = 1;
    for j = 1 : i - 1
        aux = poly_product( aux, [(1 - j) / st, 1] );
        fac = fac / (i - j) * st;
    end
    
    for j = i + 1 : st + 1
        aux = poly_product( aux, [(1 - j) / st, 1] );
        fac = fac / (i - j) * st;
    end
    
    for j = 1 : st + 1
        aux(j) = aux(j) * (1 / l) ^ (j - 1);
    end
    
    if i > 1 && i < st + 1
        pol(i + 1, :) = fac * aux;
    end
    if i == 1
        pol(i, :) = fac * aux;
    end
    if i == st + 1
        pol(2, :) = fac * aux;
    end
end

end