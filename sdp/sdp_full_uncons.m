function [SDP] = sdp_full_uncons(y, Q, x_star)
%SDP_FULL_UNCONS: create an SDP in SDPT3 format for the distance to inverse
%optimal control problem
%
%alpha_n is not eliminated in this formulation (will be in the next one)
%
% f_i(x) = (x-x_star{i})'*Q{i}*(x-x_star{i})
% f_alpha = sum_i f_i(x)
%
%INPUTS (may change around)
%   y:      test point outside the GO set
%   Q:      Dictionary of hessians
%   x_star: Corresponding optimal points to quadratic functions
%
%OUTPUTS 
%   SDP:    output SDP structure
%
%% create indexers and define variables
%form and index out the moment matrix. remember    


n = length(x_star{1});
m = length(x_star);

M_dim = 1 + n + m;

% x_index = 1 + (1:n);
% alpha_index = (n+1) + (1:(m));


%% get the blocks
blk = {'s', [M_dim]};


%% create the cost function
C = sdp_distance_cost(y, n, m);

%% Create the constraints


%% package up the SDP
% SDP = struct('blk', blk, 'C', C); %this kind of initialization forms
% struct arrays
SDP = struct;
SDP.blk = blk;
SDP.C = C;
end
