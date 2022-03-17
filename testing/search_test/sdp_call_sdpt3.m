%% only call the SDP, no rounding.
%generate the sdp

load('opt_5_2d.mat')

SOLVE = 1;

y_out = [2; 1];

[SDP_out,recoverdata,diagnostic,interfacedata] = sdp_full_uncons_yalmip(y_out, Q, x_star);

if SOLVE
[obj,X,y,Z,info,runhist] = sdpt3(SDP_out.blk,SDP_out.A,SDP_out.C,SDP_out.b); %,OPTIONS,X0,y0,Z0)

dist_rec = sqrt(obj(2));
% M = Z{end};
M = X{end};

% save('sdp_opt_5.mat', 'SDP_out', 'M', 'Z');
save('sdp_opt_5_dual.mat', 'SDP_out', 'M', 'X', 'y', 'Z', 'obj');
end