function [Erg] = energy_total(X1, KK1, midpoints_pairs, vertices2, edges2, kappa, d, gamma, points, external_force)
% returns the total energy of the problem (needed for the gradient decent algorithm)

%inputs: solution X1 of FEM system, the matrix KK1 of FEM system,
% Nx2 array of pairs of indices where midpoints (chromosomes) are located,
% array of vertices, array of edges, coefficients kappa, d, gamma,
% points at which external force is applied and the value of that force

Erg0 = 0.5 * (KK1 * X1)' * X1; %energy coming from elasticity of bundles
Erg1 = 0; %energy coming from elasticity of chromosomes
minnorm = 1e30;
for j = 1 : size(midpoints_pairs,1)
	 index1 = midpoints_pairs(j,1);
	 if index1 ~= size(vertices2,1)
		 u1 = X1(3*(index1-1)+1:3*(index1-1)+3);
	 else
		 u1 = [0,0,0]';
	 end        
	 index2 = midpoints_pairs(j,2);
	 v2 = vertices2(index2,:)';
	 v1 = vertices2(index1,:)';
	 u2 = X1(3*(index2-1)+1:3*(index2-1)+3);
	 vec = v2 - v1 + u2 - u1;
	 if norm(vec) < d
		Erg1 = Erg1 + 0.5 * kappa*(d-norm(vec))^(gamma+1); 
	 end
	 if norm(vec)<minnorm
		 minnorm = norm(vec);
	 end
end
Erg2 = 0; %energy coming from external force
for i = 1 : size(points,2)
	position = vertices2(points(i),3) + X1(3*(points(i)-1)+3);
	Erg2 = Erg2 + 0.5 * external_force * position^2;
end
Erg = Erg0 + Erg1 + Erg2;    


return