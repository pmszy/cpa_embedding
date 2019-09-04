# cpa_embedding
code for cpa(tmt) embedding

1. The current setting in module_Global.for is for the cpa_embedding case.
If you want to do tmt embedding, turn on switches mb_tmt and mb_tmt_ver1.

2. In the input file in1.data, the paramters that are relevent are
cluster size: should be L^3 for since we use simple cubic
meas, run(set to 1), so that the total number of disorder configurations is n_cores*meas*run
niter: # of iterations
nover: parameter for k-mesh
delta: energy range from -delta to delta
nwn: wn from -nwn to nwn
tables_type: set to 3 for cubic lattice
taa1: nn hopping
tz: = taa1
V1a: disorder strength
liz: linear size of liz

3. The code initialize the liz lattice by reading the file lizList, which is created by a python script lizList.py.
In lizList.py, change the value of L(linear size of the cluster) and liz and run the python script

