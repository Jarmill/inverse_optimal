function h = plot_feasible_set_constrained(A,b,C,d,xl,xu,varargin)
%PLOT_FEASIBLE_SET_CONSTRAINED plots the feasible set of a linearly
%constrained program.
%   
%   xl <= x <= xu
%   A * x == b
%   C * x <= d
%   
%   h = PLOT_FEASIBLE_SET_CONSTRAINED(A,b,C,d,xl,xu)
%   h = PLOT_FEASIBLE_SET_CONSTRAINED(A,b,C,d,xl,xu,'Name','Value')
%   
%   Inputs:
%   A,b,C,d,xu,xl ~ constraint parameters (Linear equality matrices, linear
%   inequality matrices, bound constraints vectors)
%   'Name', 'Value' ~ will be retransmitted to the plotting function
%   
%   Outputs:
%   h ~ plot handle of the feasible set

% Only bounds must be nonempty
if isempty(xl) || isempty(xu)
    error('Bound constraints must be passed.')
end
% Must be of the same length
if ~isequal(length(xl), length(xu))
    error('Bound constraints must be of same length.')
end
% Column vectors
xl = reshape(xl, [], 1);
xu = reshape(xu, [], 1);

% Determine the dimensionality of the problem
n = size(xl, 1);

% Only accept 2 or 3 dimensional problems
if n~=2 && n~=3
    error('Only 2D or 3D problems are accepted.')
end

% Input check: Check the dimensions of other inputs
if (~isempty(A) && size(A, 2)~=n) || (~isempty(b) && ~isequal(size(A, 1), size(b,1))) || (~isempty(C) && size(C, 2)~=n) || (~isempty(d) && ~isequal(size(C,1), size(d,1)))
    error('A,C must have the same number of columns as the dimension of the problem. A,C must have the same number of rows as b,d (resp.).')
end

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

% Plot differently for different dimensionality
% For 2 dimensional
if n == 2
    % If number of vertices is just 2, draw a line
    if size(V, 2) == 2
        h = plot(V(1, :), V(2, :), varargin{:});
        return
    % If number of vectices is just 1, draw a point
    elseif size(V, 2) == 1
        h = plot(V(1, :), V(2, :), 'o', varargin{:});
        return
    end
    % Else draw the convex hull of the shape
    CH = convhull(V(1, :), V(2, :));
    pgon = polyshape(V(1, CH), V(2, CH));
%     h = patch('Faces', CH.', 'Vertices', V, 'FaceColor', cool(1), 'FaceAlpha', 0, varargin{:});
    h = plot(pgon, 'FaceColor', cool(1), 'FaceAlpha', 0, varargin{:});
    return
    
% For 3 dimensional
elseif n == 3
    % If number of vertices is just 2, draw a line
    if affine_dimension(V, eps) == 1
        h = plot3(V(1, :), V(2, :), V(3, :), varargin{:});
        return
    % If number of vertices is just 1, draw a point
    elseif affine_dimension(V, eps) == 0 
        h = plot3(V(1, :), V(2, :), V(3, :), 'o', varargin{:});
        return
    % If the number of vertices is just 3, hardcode the convex hull
    elseif affine_dimension(V, eps) == 2
        % Face of the convex hull is defined just by the 3 vertices
        CH = [1 2 3 1].';
        % Plot the surface of the convex hull
%         h = trisurf(CH, V(1, :), V(2, :), V(3, :), 'FaceColor', cool(1), 'FaceAlpha', 0);
%         h = fill3(V(1, CH), V(2, CH), V(3, CH), 'FaceColor', cool(1), 'FaceAlpha', 0);
        h = fill3(V(1, CH), V(2, CH), V(3, CH), cool(1), 'FaceAlpha', 0);
        return
    end
    % Else draw the convex hull of the shape
    CH = convhull(V(1, :), V(2, :), V(3, :),'Simplify',true);
    % Remove lines as it make the view sketchy
    h = trisurf(CH, V(1, :), V(2, :), V(3, :), 'FaceColor', cool(1), 'FaceAlpha', 0, 'LineStyle', '-', varargin{:});
    return
end

end