clear all;
close all;
clc;

% Define polyhedron in 2D
n = 3;

% Define lower bounds
Alb = -eye(3);
blb = zeros(3, 1);

% Define upper bounds
Aub = eye(3);
bub = ones(3, 1);

% Define x1 + x2 <= 1
db = 1e-4;
Acons = ones(1, 3);
bcons = 1;

% Concatenate
A = [Alb; Aub; Acons];
b = [blb; bub; bcons];

% % Define x1 + x2 == 1
% db = 1e-4;
% Acons1 = ones(1, 3);
% bcons1 = 1+db;
% Acons2 = -ones(1, 3);
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
if size(V, 2) > 3
    k = convhull(V(1, :), V(2, :), V(3, :));
else
    k = [1,2,3];
end
% Figure
figure;
hold on;
hp = plot3(V(1, :), V(2, :), V(3, :), '*');
hs = trisurf(k, V(1, :), V(2, :), V(3, :), 'FaceColor', cool(1), 'FaceAlpha', 0.6);
view(3)
% patch('Faces', k.', 'Vertices', V.','FaceColor',winter(1),'FaceAlpha',.5)