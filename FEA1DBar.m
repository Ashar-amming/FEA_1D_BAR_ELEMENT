% This is a 1D FE code 
% based on uniaxial elastic two-noded bar formulation

% clear workspace
clc
clear

% read node number and coordinates
Nodes = load('Bar_Nodes.txt');
[N,l] = size(Nodes);

% read element material id, cross-sectional area
%      and connectivity
Elems = load('Bar_Elements.txt');
[E,l] = size(Elems);

% read material info
Mats = load('Bar_Materials.txt');
[M,l] = size(Mats);

% read Dirichlet and Neumann BCs
DBC = load('Bar_DBC.txt');
[P,l] = size(DBC);

NBC = load('Bar_NBC.txt');
[Q,l] = size(NBC);

% determine total number of degrees-of-freedom
udof = 1;       % dof per node
NDOF = N*udof;  % total dof for problem

% initialize matrices and vectors
K = zeros(NDOF,NDOF);   % global stiffness matrix
U = zeros(NDOF,1);      % global displacement vector
F = zeros(NDOF,1);      % global external force vector

% loop of bar elements
for e=1:E
    
    % establish element connectivity and coordinates
    Nnums = Elems(e,4:5);
    x = Nodes(Nnums(:),2);
    
    % extract element cross-sectional area
    A = Elems(e,3);
    
    % extract element Young's modulus
    Y = Mats(Elems(e,2),2);
    
    % determine element length
    L = sqrt( (x(2)-x(1))*(x(2)-x(1)) );
    
    % construct element stiffness matrix
    k = Y*A/L;
    Ke = k*[ 1 -1; -1 1 ];
    
    % assemble element stiffness into global matrix
    c1 = Nnums(1);
    c2 = Nnums(2);
    K(c1,c1) = K(c1,c1) + Ke(1,1);
    K(c1,c2) = K(c1,c2) + Ke(1,2);
    K(c2,c1) = K(c2,c1) + Ke(2,1);
    K(c2,c2) = K(c2,c2) + Ke(2,2);
    
end

% construct global force vector
for q=1:Q
    F(NBC(q,2)) = NBC(q,3);
end

% set penalty for displcement constraints
Klarge = 10^6;

% enforce displacement constraints
for p=1:P
    K(DBC(p,1),DBC(p,1)) = Klarge;
    F(DBC(p,1)) = Klarge*DBC(p,2);
end

format long
K
F

% solve system for nodal displacements
U = K\F;

% output displacements
U

% recover element internal results

Disp = zeros(E,3);
Eps  = zeros(E,1);
Sig  = zeros(E,1);

for e=1:E
    
    % establish element connectivity and coordinates
    Nnums = Elems(e,4:5);
    x = Nodes(Nnums(:),2);
    
    % extract element cross-sectional area
    A = Elems(e,3);
    
    % extract element Young's modulus
    Y = Mats(Elems(e,2),2);
    
    % determine element length
    L = sqrt( (x(2)-x(1))*(x(2)-x(1)) );
    
    % extract element nodal displacements
    Disp(e,1) = U(Nnums(1));
    Disp(e,2) = U(Nnums(2));
    Disp(e,3) = ( U(Nnums(1))+U(Nnums(2)) )/2;
    
    % calculate element strain
    Eps(e,1) = ( Disp(e,2) - Disp(e,1) )/L;
    
    % calculate element stress
    Sig(e,1) = Y*Eps(e,1);
    
end

Disp
Eps
Sig
