This folder may be unnecessary.

YALMIP does get the moment information, stored in the dual variable Z
associated with the PSD matrix (Z{end}).

There is a naming conflict between YALMIP and SDPNAL in the term 'constraint'

I'm not sure how to address this. It may be best to precompile the SDPs in YALMIP
for each aspect of the dataset, and then turn give the SDPs to STRIDE.