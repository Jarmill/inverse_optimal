function V = get_vertices_of_polytope(A, b, C, d, xl, xu)
%GET_VERTICES_OF_POLYTOPE finds the vertices of a polytope described by
%linear constraints.

% Find the number of variables
n = length(xl);

% Define a slackness for the equality constraints
db = 1e-6;
% Create the polyhedron representation: Apoly * x <= bpoly
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

% Merge them all together into polyhedron inequality constraints
Apoly = [Alb; Aub; Aec1; Aec2; Aic;];
bpoly = [blb; bub; bec1; bec2; bic;];

% Get the vertices of this polyhedron
[V, ~] = con2vert(Apoly, bpoly);
% Get the vertices of this polyhedron into column vector format
V = V.';

% Since the equality constraints are given slack, we'll have additional
% vertices. Remove those closer than 2*sqrt(n)*db from each other.
V = reduce_vertices(V, 2*sqrt(n)*db);


end