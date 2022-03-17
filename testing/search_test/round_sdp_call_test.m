load('opt_5_2d.mat')

y_out = [2; 1];

[SDP_out,recoverdata,diagnostic,interfacedata] = sdp_full_uncons_yalmip(y_out, Q, x_star);

[obj,X,y,Z,info,runhist] = sdpt3(SDP_out.blk,SDP_out.A,SDP_out.C,SDP_out.b); %,OPTIONS,X0,y0,Z0)

dist_rec = sqrt(obj);
M = Z{end};

save('sdp_opt_5.mat', 'SDP_out', 'M', 'Z');