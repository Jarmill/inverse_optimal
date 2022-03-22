function x = global_optima_set_qp(Q, phi, A, b, C, d, xl, xu, Nmesh)
%GLOBAL_OPTIMA_SET_QP finds a regular meshgrid of points describing the set
%of global optima, for a general convex QP.
%
%   min_x    1/2*x'*Q*x + phi'*x
%   s.t.     xl <= x <= xu
%            A*x == b
%            C*x <= d
%   
%   Notes: Equality constraints are decomposed as inequality constraints
%   A*x <= b+db& -A*x <= -(b-db)
%   with small db.
%   
%   x = global_optima_set_qp(Q, phi, A, b, C, d, xl, xu, Nmesh)
%   
%   Inputs:
%   Q, phi, A, b, C, d, xl, xu ~ Quadratic Program parameters
%   Nmesh ~ maximum number of mesh points to generate from the polyhedron
%   of solutions to the QP.
%   
%   Outputs:
%   x ~ matrix of grid-points from the polyhedron of solutions to the QP

% Extract size of QP
n = size(Q,1);

% Consistency of inputs
if ~isequal(size(Q,2), n) || ~isequal(size(phi, 1), n) || (~isequal(A, []) && ~isequal(size(A, 2), n)) || ...
    (~isequal(C, []) && ~isequal(size(C, 2), n)) || ~isequal(size(xl, 1), n) || ~isequal(size(xu, 1), n)
    error('Phi, xl, xu should be column vectors of size equal to the variable of the QP. Q must be square. Q, A, C should have the same number of columns as size of variable of the QP.');
end
if ~isequal(size(A, 1), size(b, 1)) || ~isequal(size(C, 1), size(d, 1))
    error('A and b, C and d, should have the same number of rows.');
end

% Solve the QP and get a single solution
qpopt = optimoptions(@quadprog, 'Algorithm', 'Interior-point-convex', 'Display', 'Iter');
[x_opt, ~, exitflag, ~] = quadprog(Q, phi, C, d, A, b, xl, xu, [], qpopt);

% If anything but convergence, or step-tolerance issue a warning
if ~isequal(exitflag, 1) && ~isequal(exitflag, 2)
    warning('global_optima_set_qp: The solution that was found is not guaranteed to be an optimum, thus the mesh output by the function may be erroneous.');
end    

% Now that you have a single optimum of the QP you can formulate the
% polyhedron representation of its set of solutions

% If Q is full rank don't do that though, return the optimal solution
if rank(Q) == n
    x = x_opt;
    return
end

% If Q is not full rank proceed to the polyhedron formulation

% Get orthonormal basis for range of Q
RQ = orth(Q);
% Get orthonormal basis for nullspace of Q
NQ = null(Q);

% Project phi onto the range and nullspace
phiR = RQ * RQ.' * phi;
phiN = NQ * NQ.' * phi;

% LP formulation of a convex QP is
%
% min_x     phiN'*x
% s.t.      xl <= x <= xu
%           A*x == b
%           C*x <= d
%           Q*x == -phiR
%
% Thus if we have an optimal solution x_opt, we can characterize all the
% solutions of the LP as a polyhedron
% 
% x in S <=> x in R^n
% s.t.      xl <= x <= xu
%           A*x == b
%           C*x <= d
%           Q*x == -phiR
%           phiN'*x <= phiN'*x_opt
%
% Characterize the polyhedron in terms of inequality constraints
%   
%   Notes: Equality constraints are decomposed as inequality constraints
%   A*x <= b+db & -A*x <= -(b-db)
%   Q*x <= -phiR+db & -Q*x <= -(-phiR-db)

% Define the equality constraint slackness db
db = 1e-5;

% Lower and upper bound constraints
Alb = -eye(n);
blb = -xl;
Aub = eye(n);
bub = xu;
% Equality constraints are a pair of inequalities
Aec1 = A;
bec1 = b+db;
Aec2 = -A;
bec2 = -(b-db);
% Inequality constraints
Aic = C;
bic = d;
% Partial optimality conditions are equality constraints, so a pair of
% inequalities
Aoc1 = Q;
% boc1 = -phiR+db;
boc1 = Q*x_opt+db;
Aoc2 = -Q;
% boc2 = -(-phiR-db);
boc2 = -(Q*x_opt-db);
% Minimum conditions are inequality constraints 
Amc = phiN.';
bmc = phiN.'*x_opt+db;

% Merge them all together into polyhedron inequality constraints
Apoly = [Alb; Aub; Aec1; Aec2; Aic; Aoc1; Aoc2; Amc];
bpoly = [blb; bub; bec1; bec2; bic; boc1; boc2; bmc];

% % Check feasibility
% num_conspoly = [length(blb), length(bub), length(bec1), length(bec2), length(bic), length(boc1), length(boc2), length(bmc)];
% name_conspoly = ["lb", "ub", "ec1", "ec2", "ic", "oc1", "oc2", "bmc"];
% 
% cumnum_conspoly = cumsum(num_conspoly)
% 
% val_conspoly = Apoly*x_opt-bpoly;
% ind_infeasible = find(val_conspoly >= 0)
% infeasibleval = val_conspoly(ind_infeasible)
% for infinf = 1 : length(ind_infeasible)
%     fprintf('%s\n', name_conspoly(find(cumnum_conspoly >= ind_infeasible(infinf), 1, 'first')));
% end
% 
% % Is solvable
% lpopt = optimoptions(@linprog, 'Display', 'Iter');
% [~, ~, exitflaglp, ~] = linprog(zeros(n, 1), Apoly, bpoly, [], [], [], [], lpopt);


% Get the vertices of this polyhedron
[V, nr] = con2vert(Apoly, bpoly);
% Get the vertices of this polyhedron into column vector format
V = V.';

% Since the equality constraints are given slack, we'll have additional
% vertices. Remove those closer than 2*sqrt(n)*db from each other.
V = reduce_vertices(V, 2*sqrt(n)*db);

%In order to produce a regular grid over the polyhedron, generate a regular
%grid of convex combinations of the polyhedron vertices.
%Convex combination parameters lie on a simplex. So generate a grid on the 
%simplex.
alpha = prob_simplex_ndim(size(V, 2)-1, Nmesh);
% Get the parametrizations into column vector format
alpha = alpha.';

% Get the grid of solutions to the QP as a grid over the polyhedron of
% solutions.
x = V * alpha;
end