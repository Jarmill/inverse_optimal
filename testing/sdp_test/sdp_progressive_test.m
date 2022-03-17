%gradually build up the SDP for uncons optimization

load('opt_5_2d.mat')

y_out = [2; 1];

SDP = sdp_full_uncons(y_out, Q, x_star);

% [obj,X,y,Z,info,runhist] = sdpt3(SDP_out.blk,SDP_out.A,SDP_out.C,SDP_out.b); %,OPTIONS,X0,y0,Z0)