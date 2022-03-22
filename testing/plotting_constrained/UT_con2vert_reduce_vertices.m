clear all;
close all;
clc;

% Define polyhedron in 2D
n = 2;

% Define lower bounds
Alb = -eye(2);
blb = zeros(2, 1);

% Define upper bounds
Aub = eye(2);
bub = ones(2, 1);

% Define x1 + x2 <= 1
db = 0.0001;
Acons1 = ones(1, 2);
bcons1 = 1;

% Concatenate
A = [Alb; Aub; Acons1];
b = [blb; bub; bcons1];

% % Define x1 + x2 == 1
% db = 0.0001;
% Acons1 = ones(1, 2);
% bcons1 = 1+db;
% Acons2 = -ones(1, 2);
% bcons2 = -(1-db);
% 
% % Concatenate
% A = [Alb; Aub; Acons1; Acons2];
% b = [blb; bub; bcons1; bcons2];

% Find the vertices
[V, nr] = con2vert(A, b);
% Transform vertices in column vector format
V = V.';
% Reduce the vertices
V = reduce_vertices(V, 1.2*db);
% Find the number of vertices
nv = size(V, 2);
% Find the convex hull of the vertices
if size(V, 2) > 2
    k = convhull(V(1, :), V(2, :));
else
    k = [1;2];
end
% Figure
figure;
hold on;
plot(V(1, :), V(2, :), '*')
patch('Faces', k.', 'Vertices', V.','FaceColor',winter(1),'FaceAlpha',.5)