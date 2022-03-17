%% Example: solve the ProjGO SDP relaxation using STRIDE
%this script is based on the QUASAR example from STRIDE

%% Start clean
% clc; clear; close all; restoredefaultpath;

%% Add STRIDE to matlab path, and provide path to dependencies
% addpath(genpath(pwd));
manoptpath      = '../manopt'; % required for local search
sdpnalpath      = '../SDPNAL+v1.0'; % required for ADMM+

%% Load problem data

load("sdp_opt_5.mat")
load('opt_5_2d.mat')


y = [2; 1];


%% Solve using STRIDE
addpath(genpath(manoptpath)); % add manopt to path

% set parameters for STRIDE
options.pgdStepSize     = 10; % step size, default 10
options.maxiterPGD      = 5; % maximum outer iterations for STRIDE, default 5-10
options.SDPNALpath      = sdpnalpath; % provide path to SDPNAL
options.tolADMM         = 1e-4; % tolerance for warmstart, decrease this parameter for a better warmstart (but takes more time)
options.tolPGD          = 1e-8; % tolerance on KKT residual of the SDP
options.lbfgseps        = false;


% provide implementation to the local search method
options.rrOpt           = 1:3; % round the leading 3 eigenvectors to generate hypotheses
options.rrFunName       = 'local_search_uncons'; % name of the .m file that implements the local search
% options.rrPar = 



rrPar = struct;
rrPar.Q = Q;
rrPar.x_star = x_star;
rrPar.y = y;
% rrPar.m = m;
% Primal initialization
X0                  = [];

% options.printlevel      = false;

% call STRIDE
[outPGD,Xopt,yopt,Sopt] = PGDSDP(SDP_out.blk, SDP_out.A, SDP_out.b, SDP_out.C, X0, options);
infostride              = get_performance_quasar(Xopt,yopt,Sopt,SDP,R_gt);



