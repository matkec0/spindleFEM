function [S] = FEM_matrix(vertices, edges, QHQT, st1, st2, penal)
% FEM matrix for system of elastic rods given for graph given with 
% array of vertices and array of edges;
% QHQT is array of matrices, equal to product of matrices Q*H*Q^T, 
% important for the model of the elastic rod; st1 and st2 are desired
% powers (k in Pk elements) for discretization FEM space; penal
% is probably not needed



global nb % array of normal vectors

vrN = size(vertices,1);
edN = size(edges,1);


AFORMA = a_forma(st1, 1); %% u'*v', za 1. clan
B1FORMA = b1_forma(st1, st2, 1); %% u' * v¸, za 4. i 5. clan
B2FORMA = b2_forma(st1, st2, 1); %% u * v, za 4. i 5. clan
CFORMA = fixing(st1, 1);
CFORMA2 = kron(CFORMA,eye(3)); %%% za fiksiranje

%indices where first occurence of displacement u, infitesimal rotation om, 
%	Lagrange multiplier lambda, and fixing parameters a and b, respectively
%	are located in the solution vector
first_u = 1;
first_om = vrN * 3 + edN * (st1 - 1) * 3 + 1;
first_lambda = vrN * 6 + edN * (st1 - 1) * 6 + 1;
first_aib = vrN * 6 + edN * (st1 - 1) * 6 + edN * (st2+1)*3 + 1;

% dimension of the system
m = vrN * 6 + edN * ( (st1 - 1) * 6 + (st2 + 1) * 3 ) + 6;
S = sparse(m,m);


for ii = 1 : edN
    % taking all important parameters for building one elementary
	% matrix, for one edge
	edge = edges(ii, :);
	v1 = vertices(edge(1), :);
	v2 = vertices(edge(2), :);
	l = norm(v2 - v1);
    t = (v2 - v1)' / l; 
    LQHQT = QHQT(3*(ii-1)+1:3*(ii-1)+3,1:3) / l; % 1/L * QHQT 
	At = -[0, -t(3), t(2); t(3), 0, -t(1); -t(2), t(1), 0]; 
	LAt = l*At;
    
	
	eye3=eye(3);
	E=zeros(6*(st1+1)+3*(st2+1)); % one elementary matrix for the edge
	E(3*(st1+1)+1:6*(st1+1),3*(st1+1)+1:6*(st1+1)) = kron (AFORMA,LQHQT);
	E(1:3*(st1+1),6*(st1+1)+1:6*(st1+1)+3*(st2+1))=	kron (B1FORMA,eye3);
	E(3*(st1+1)+1:6*(st1+1),6*(st1+1)+1:6*(st1+1)+3*(st2+1)) = kron (B2FORMA,LAt);
	E(6*(st1+1)+1:6*(st1+1)+3*(st2+1),1:3*(st1+1))=	kron (B1FORMA,eye3)';
	E(6*(st1+1)+1:6*(st1+1)+3*(st2+1),3*(st1+1)+1:6*(st1+1)) = kron (B2FORMA,LAt)';


	% now we have to find where to put the elementary matrix E in the whole matrix S
	% for this we need the array "indices", denoting those positions
	indices = zeros (2*(st1+1)+(st2+1),1);
	indices(1)= first_u  + (edge(1) - 1) * 3;
	indices(2)= first_u  + (edge(2) - 1) * 3;
	for j=3:(st1+1) %the rest of indices for u
		indices(j) = first_u + vrN * 3 + (ii - 1) * 3 * (st1 - 1) + 3*(j-3);
	end
	indices((st1+1)+1)= first_om  + (edge(1) - 1) * 3;
	indices((st1+1)+2)= first_om  + (edge(2) - 1) * 3;
	for j=3:(st1+1) %the rest of indices for om
		indices((st1+1)+j) = first_om + vrN * 3 + (ii - 1) * 3 * (st1 - 1) + 3*(j-3);
	end	
	for j=1:(st2+1) %all indices for lambdas
		indices(j+2*(st1+1)) = first_lambda + (ii - 1) * 3 * (st2 + 1) + 3*(j-1);
	end

	% now we transfer E to S
	for j1=1:length(indices)
		for j2=1:length(indices)
			tempM = E((3*j1-2):(3*j1),(3*j2-2):(3*j2));
			if norm(tempM,'fro')~=0
				S(indices(j1):indices(j1)+2,indices(j2):indices(j2)+2) = S(indices(j1):indices(j1)+2,indices(j2):indices(j2)+2)+tempM;
			end
		end
	end

	% positions regarding the unknowns a and b need to go to last six rows and columns of the matrix S
	for j1=1:(st1+1)
		tempM=l*CFORMA2*penal;
		S(indices(j1):indices(j1)+2,first_aib:first_aib+2) = S(indices(j1):indices(j1)+2,first_aib:first_aib+2) + tempM(3*j1-2:3*j1,:);
		S(first_aib:first_aib+2,indices(j1):indices(j1)+2) = S(first_aib:first_aib+2,indices(j1):indices(j1)+2) + tempM(3*j1-2:3*j1,:)';
	end
	for j1=1:(st1+1)
		j2=j1+(st1+1);
		tempM=l*CFORMA2*penal;
		S(indices(j2):indices(j2)+2,first_aib+3:first_aib+5) = S(indices(j2):indices(j2)+2,first_aib+3:first_aib+5) + tempM(3*j1-2:3*j1,:);
		S(first_aib+3:first_aib+5,indices(j2):indices(j2)+2) = S(first_aib+3:first_aib+5,indices(j2):indices(j2)+2) + tempM(3*j1-2:3*j1,:)';
	end
end

end