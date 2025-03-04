function [r,theta] = r_theta(N)
% returns pairs of polar coordinates for N shells

vec_a = zeros (N*(N+1)/2,1);
vec_b = zeros (N*(N+1)/2,1);
count = 0;
for i = 1 : N
	for j = 1 : i
		count = count +1;
		vec_a(count) = j;
		vec_b(count) = i-j;
	end
end


[r,theta] = r_theta_a_b (vec_a, vec_b);
X = [r,theta];
X = sortrows(X);
r = X(:,1);
theta = X(:,2);
end

function [r,theta] = r_theta_a_b (a, b)
% for a point P that that is on the position a*e1 + b *v
% (v is a vector in hexagon latice pointing to up-right), returns
% r ... distance of P to the origin
% theta ... polar angle
% if a and b are vectors, then return value is Nx2 matrix


N = length(a);
r 		= zeros(N,1);
theta 	= zeros(N,1);
for i = 1 : N
	r(i) 	 = sqrt (a(i)^2 + b(i)^2 + a(i)*b(i));
	theta(i) = acos ((2*a(i)+b(i))/(2*r(i)));
end 

end