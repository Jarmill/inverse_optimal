function [model,indexer, info] = uncons_export(y, Q, x_star, solver)
%SDP_FULL_UNCONS: create an SDP in SDPT3 format for the distance to inverse
%optimal control problem
%
% f_i(x) = (x-x_star{i})'*Q{i}*(x-x_star{i})
% f_alpha = sum_i f_i(x)
%
%INPUTS (may change around)
%   y:      test point outside the GO set
%   Q:      Dictionary of hessians
%   x_star: Corresponding optimal points to quadratic functions
%
%Outputs 
%   model:  Exported SDP model
%   indexer:index and output structure for yalmip
%   info:   details of the YALMIP export

if nargin <= 4
    solver = 'sdpt3';
end
opts = sdpsettings('solver', solver);
    

[Fd, objd, indexer] = sdp_full_uncons_yalmip(y, Q, x_star);

[model,recoverdata,diagnostic,interfacedata] = export(Fd, -objd, opts);

info = struct;
info.recoverdata = recoverdata;
info.diagnostic= diagnostic;
info.interfacedata= interfacedata;
end