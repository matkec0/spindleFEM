%%%% main MATLAB script for numerical scheme modelling the stucture 
%%%% of spindle and chromosomes
%%%% output is array new_positions, denoting new positions of
%%%% spindle points


close all
clear
addpath('./auxiliaryFunctions/');
format shortE


% auxiliary variables for the case of clamped boundary conditions in poles
clamped = 0;
anglex = 70;
maxangle = anglex*pi/180;
omega3 = tan(maxangle);


N_shells = 6; %number of shells, each shells containing N_bundles chromosomes
N_bundles = 6; %number of chromosomes in each shells
N_spindle = N_shells * N_bundles; %total number of bundles (up to 1, central one, will be incorporated later)
% total number of chromosomes is N_shells*N_bundles+1 (1 is on the central line, will be incorporated later)
L = 14;   %spindle length
d = 3; %chromosome diameter
nuCHR = 0.3;
ECHR = 1000; % elastic coefficients for chromosomes
CR = @(radijus) 1/(1.5 * (1-nuCHR^2)/ECHR*sqrt(2/radijus));
C0 = CR(d/2.);  
kappa = 4*C0/5; %constants in the definition of chromosome repulsive forces
gamma = 1.5; %the power constant in chromosome repulsive force

kappa0 = 1e-1;  %numerics parameter: the greater, the convergence is slower, but more stable
N_FEM = 14; %number of discretization points on each spindle
MP_shift = 0; % metaphase plane shift (in number of discretization points); if =0, the plane is in the middle
external_force = 0;  %the value of external force applied to the whole spindle

tol = 1e-8; %tolerance for breaking the gradient descent loop
N_steps = 100000; %maximal number of steps in gradient descent
alpha = min(1e-1,5e1/C0^2); %step length in gradient descent


%%% building bundles and other auxiliary variables

edges0 = [1:N_FEM;2:N_FEM+1]'; % auxiliary edges in one bundle
vertices2 = [];
edges2 = [];
nb0 = zeros(N_FEM,3);
nb0(:,2)=ones(N_FEM,1); %normal vectors
nb = []; %normal vectors
Q = []; % Frenet frames

[r,thetas0] = r_theta (5); % ako je u zagradi N, onda imamo ukupno 
% r ... distances of first chromosomes in each shell in metaphase plane to the pole (projected to the plane)
% thetas0 ... angles of first chromosomes in each shell in metaphase plane to the pole (projected to the plane)
% for N=5 argument we have in total N*(N+1) * 3 "+ 1" chromosomes; can be more than N=5 

if clamped == 1
    angle0 = (r/r(N_shells))*omega3;
    omeg0 = [];
end

for i = 1 : N_shells 
    vL = r(i)*[0,1,0] + [-L/2,0,0];
    vR = r(i)*[0,1,0] + [ L/2,0,0];
    distance = vR-vL;
    curve0 = [vL(1):distance(1)/N_FEM:vR(1);r(i)*ones(1,N_FEM+1);zeros(1,N_FEM+1)]';
    bn0 = zeros(N_FEM,3);
    bn0(:,3)=ones(N_FEM,1);
    t0 = zeros(N_FEM,3); %binormals and tangents, for Frenet frames
    Q0=[];
    for j = 1 : N_FEM
        t0(j,:) = (curve0(j+1,:)-curve0(j,:))/norm(curve0(j+1,:)-curve0(j,:));
        nb0(j,:) = cross(bn0(j,:)',t0(j,:)')';
        Q0 = [Q0;[t0(j,:);nb0(j,:);bn0(j,:)]];
    end
	% vertices (curve0), edges and Frenet frames for one rod in shell created 
	for_angle = thetas0(i) : 2*pi/N_bundles : 2*pi;
    for ii = 1 : N_bundles
        rot = [[1,0,0];[0,cos(for_angle(ii)),sin(for_angle(ii))];[0,-sin(for_angle(ii)),cos(for_angle(ii))]];
        spindlex = (rot*curve0')';
        vertices2 = [vertices2;spindlex];
        edges2=[edges2;edges0+(N_FEM+1)*(ii-1)*ones(size(edges0))+(i-1)*(N_FEM+1)*N_bundles*ones(size(edges0))];
        nbx = (rot*(nb0'))';
        nb = [nb;nbx];
        Qx = (rot*Q0')';
        Q = [Q;Qx];
        if clamped == 1
            omeg0 = [omeg0;(rot*[0;0;angle0(i)])'];
        end
    end
end
if clamped == 1
    omeg0 = [omeg0;-omeg0];
end
%%% end of building bundles and other auxiliary variables

vrN = size(vertices2,1);
edN = size(edges2,1);
shift_rot = 3*vrN+3*edN; % shift of infitesimal rotation in solution of FEM system

midpoints = [vrN+1, N_FEM/2+1 : N_FEM+1 : vrN-1]; %indices of midpoints of bundles
midpoints = midpoints + MP_shift*ones(size(midpoints)); %shifting in case of nontrivial MP_shift
midpoints(1)=vrN+1;
midpoints_pairs = []; %pairs of all midpoints (if too close, chromosome repulsive forces are activated)
for i = 1 : length(midpoints)-1
    for j = i+1 : length(midpoints) 
        midpoints_pairs = [midpoints_pairs; midpoints(i), midpoints(j)];
    end
end

vr2plus=[vertices2;vertices2(midpoints(2),1),0,0];
% adding the last vertex (midpoint of the central line)


% QHQT are important matrices in elastic rod model (for bundles)
% here those matrices are formed
H = zeros(3, 3);
QHQT = zeros (3*edN,3);
H(1, 1) = 900;
H(2, 2) = 900;
H(3, 3) = 900;
for ii = 1 : edN
    indexb = 3*(ii-1);
    indexv1 = 3*(edges2(ii,1)-1);
    indexv2 = 3*(edges2(ii,2)-1);
    auxQ = Q(indexb+1:indexb+3,1:3)';
    QHQT(indexb+1:indexb+3,1:3) = 0.5 * ((auxQ * H * auxQ') + (auxQ * H * auxQ')');
end


%%% forming FEM matrix
KK1 = FEM_matrix(vertices2, edges2, QHQT, 2, 1, 100);


%% applying Dirichlet boundary conditions on displacements in poles
springs_u = [1:N_FEM+1: vrN, N_FEM+1:N_FEM+1:vrN];
II3 = eye(3);
springs_u_ = springs_u(1:N_spindle);
Kpenal = 1e20;
for ii = 1 : length(springs_u_)
    indices_i = ( 3*(springs_u_(ii)-1) + 1 ) : ( 3*(springs_u_(ii)-1) + 3 );
    KK1(indices_i, indices_i)= KK1( indices_i, indices_i ) + Kpenal*II3;
end
springs_u_ = springs_u(N_spindle+1:2*N_spindle);
for ii = 1 : length(springs_u_)
    indices_i = ( 3*(springs_u_(ii)-1) + 1 ) : ( 3*(springs_u_(ii)-1) + 3 );
    KK1(indices_i, indices_i)= KK1( indices_i, indices_i ) + Kpenal*II3;
end



%springs between midpoints of bundles
KK1k = KK1;
for j = 1 : size(midpoints_pairs,1)
    index1 = midpoints_pairs(j,1);
    index2 = midpoints_pairs(j,2);
    if index1 ~= vrN+1
        KK1k(3*(index1-1)+1:3*(index1-1)+3, 3*(index1-1)+1:3*(index1-1)+3) = KK1k(3*(index1-1)+1:3*(index1-1)+3, 3*(index1-1)+1:3*(index1-1)+3) + kappa0*eye(3);
        KK1k(3*(index1-1)+1:3*(index1-1)+3, 3*(index2-1)+1:3*(index2-1)+3) = KK1k(3*(index1-1)+1:3*(index1-1)+3, 3*(index2-1)+1:3*(index2-1)+3) - kappa0*eye(3);
        KK1k(3*(index2-1)+1:3*(index2-1)+3, 3*(index1-1)+1:3*(index1-1)+3) = KK1k(3*(index2-1)+1:3*(index2-1)+3, 3*(index1-1)+1:3*(index1-1)+3) - kappa0*eye(3);
    end
    KK1k(3*(index2-1)+1:3*(index2-1)+3,3*(index2-1)+1:3*(index2-1)+3) = KK1k(3*(index2-1)+1:3*(index2-1)+3,3*(index2-1)+1:3*(index2-1)+3) + kappa0*eye(3);
end


if clamped == 1
	%% applying Dirichlet boundary conditions on rotations in poles
    springs_u = [1:N_FEM+1: vrN, N_FEM+1:N_FEM+1:vrN];
    II3 = eye(3);
    springs_u_ = springs_u(1:N_spindle);
    Kpenal = 1e20;
    for ii = 1 : length(springs_u_)
        indices_i = (shift_rot + 3*(springs_u_(ii)-1) + 1 ) : (shift_rot + 3*(springs_u_(ii)-1) + 3 );
        KK1k(indices_i, indices_i)= KK1k( indices_i, indices_i ) + Kpenal*II3;
    end
    springs_u_ = springs_u(N_spindle+1:2*N_spindle);
    for ii = 1 : length(springs_u_)
        indices_i = (shift_rot + 3*(springs_u_(ii)-1) + 1 ) : (shift_rot + 3*(springs_u_(ii)-1) + 3 );
        KK1k(indices_i, indices_i)= KK1k( indices_i, indices_i ) + Kpenal*II3;
    end
end


NN = size(KK1,1);
F0 = zeros(NN,1);

%% setting the initial configuration, for the first iteration of gradient descent
normals = [];
for ii = 1 : edN
    indexb = 3*(ii-1);
    indexv1 = 3*(edges2(ii,1)-1);
    indexv2 = 3*(edges2(ii,2)-1);
    auxQ = Q(indexb+1:indexb+3,1:3)';
    normal = auxQ(:,2);
    normals = [normals;normal];
end

if clamped == 1
    for ii = 1 : length(springs_u)
        indices_i = (shift_rot + 3*(springs_u(ii)-1) + 1 ) : (shift_rot + 3*(springs_u(ii)-1) + 3 );
        F0(indices_i)= F0(indices_i) + Kpenal*(omeg0(ii,:)');
    end
end

%solving first FEM system
KK1k = KK1k + speye(size(KK1k)) * 1e-10; %ensuring regularity of the matrix
X1 = KK1k \ F0;
X = solutionsTransform (X1,size(vertices2,1),size(edges2,1),2);

 
%external force
fun = @(x)  - external_force *  x(3) * [0;0;1]; 
F1 = FEM_RHS_at_point(fun, vertices2, midpoints(2:size(midpoints,2)), size(KK1,1));


%calculating the initial energy of system
%gradient descent tries to decrease the energy
tildeX1 = X1;
Erg0 = energy_total(tildeX1, KK1, midpoints_pairs, vr2plus, edges2, kappa, d, gamma, midpoints(2:size(midpoints,2)), external_force);


minvec=0;
vec=d;
i = 0;
Erg_err = 1e30;
max_displacement = 1e30;
new_positions = deformation(vertices2, X); % variable that keeps the displacement of the spindle at the end of this script

%% start of the main gradient descent loop 
while i < N_steps & alpha > 1e-10 & abs(Erg_err/Erg0) > tol & max_displacement/L > tol
    i = i + 1;
	F0 = FEM_RHS_at_point(fun, new_positions, midpoints(2:size(midpoints,2)), size(KK1,1));
    % applying repulsive chromosome forces
	if mod(i,1) == 0
         minvec = 1e30;
         for j = 1 : size(midpoints_pairs,1)
             index1 = midpoints_pairs(j,1);
             v1 = vr2plus(index1,:)';
             if index1 ~= vrN+1
                 u1 = tildeX1(3*(index1-1)+1:3*(index1-1)+3);
             else
                 u1 = [0,0,0]';
             end
             index2 = midpoints_pairs(j,2);
             v2 = vr2plus(index2,:)';
             u2 = tildeX1(3*(index2-1)+1:3*(index2-1)+3);
             vec = v2 - v1 + u2 - u1;
             if norm(vec) < minvec
                 minvec = norm(vec);
             end
             if norm(vec) < d %% the force is applied only if the chromosomes are touching
                 if index1 ~= vrN+1
                     F0(3*(index1-1)+1:3*(index1-1)+3) = F0(3*(index1-1)+1:3*(index1-1)+3) - 0.5*kappa*(gamma+1)*(d-norm(vec))^gamma/norm(vec)*vec;
                 end
                 F0(3*(index2-1)+1:3*(index2-1)+3) = F0(3*(index2-1)+1:3*(index2-1)+3) + 0.5*kappa*(gamma+1)*(d-norm(vec))^gamma/norm(vec)*vec;
             end
         end
     end
	 
     forcing_old = KK1 * tildeX1;
     FF = F0 - forcing_old;
     if clamped == 1
        for ii = 1 : length(springs_u)
            indices_i = (shift_rot + 3*(springs_u(ii)-1) + 1 ) : (shift_rot + 3*(springs_u(ii)-1) + 3 );
            FF(indices_i)= [0;0;0];
        end
     end
	 % solving a FEM system in gradient descent loop
	 
     X1 = KK1k \ FF;
     increment = X1;
	 
	 %searcing for maximal gradient descent step size such that the energy decreases
     while energy_total(tildeX1 + alpha * increment, KK1, midpoints_pairs, vr2plus, edges2, kappa, d, gamma, midpoints(2:size(midpoints,2)), external_force) > Erg0
         alpha = alpha/2;
     end

     % incrementing the solution    
     tildeX1 = tildeX1 + alpha * increment;
     Ergnew = energy_total(tildeX1, KK1, midpoints_pairs, vr2plus, edges2, kappa, d, gamma, midpoints(2:size(midpoints,2)), external_force);
     Erg_err = Erg0 - Ergnew;
     max_displacement = max(abs(alpha*increment(1:shift_rot)));
     Erg0 = Ergnew;
    
     if alpha < 1e-2  
         alpha = alpha*1.1;
     end
     X = solutionsTransform (tildeX1,size(vertices2,1),size(edges2,1),2);
     new_positions = deformation(vertices2, X); 
 end

new_positions = deformation(vertices2, X); 
springs_u = [1:N_FEM+1: vrN, N_FEM+1:N_FEM+1:vrN]; 
vec = [];
omeg = [];
for ii = 1 : length(springs_u)
    indexi_0 = ( shift_rot + 3*(springs_u(ii)-1) + 1 ) : ( shift_rot + 3*(springs_u(ii)-1) + 3 );
    omegpom = tildeX1(indexi_0);
    omeg = [omeg; omegpom'];
end
 
 
 
%ploting the obtained solution

 figure(30)
 clf(30)
 hold on
 hold on
 [x,y,z]=sphere;
 x1= 0.5*d* x;
 y1= 0.5*d* y;
 z1= 0.5*d* z;
 daspect([1 1 1])
 new_positions = deformation(vertices2, X);
for i = 1 : size(edges2,1)
    edge = edges2(i, :);
    prvi = new_positions(edge(1), :);
    drugi = new_positions(edge(2), :);
    x = [prvi(1); drugi(1)];
    y = [prvi(2); drugi(2)];
    z = [prvi(3); drugi(3)];
    plot3(x, y, z, 'r-','LineWidth',2);
    hold on
end
 new_positions(vrN+1,:)=[vertices2(midpoints(2),1),0,0];		 					 
for iii = 1 : size(midpoints,2)
	surf(x1+new_positions(midpoints(iii),1),y1+new_positions(midpoints(iii),2),z1+new_positions(midpoints(iii),3),'FaceAlpha',.3,'EdgeColor','none');
end
view(90,0)
daspect([1 1 1])
pause(1)
view(-40,40)
daspect([1 1 1])
pause(1)
view(0,0)
daspect([1 1 1])
pause(1)
hold off


if clamped == 0
    omeg0 = omeg;
end

