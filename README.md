# cpa_embedding
code for cpa(tmt) embedding

1. The current setting in module_Global.for is for the cpa_embedding case.
If you want to do tmt embedding, turn on switches mb_tmt and mb_tmt_ver1.

2. In the input file in1.data, the paramters that are relevent are
1) luster size: should be L^3 for since we use simple cubic
2) meas, run(set to 1), so that the total number of disorder configurations is n_cores*meas*run
niter: # of iterations
3) nover: parameter for k-mesh
4) delta: energy range from -delta to delta
5) nwn: wn from -nwn to nwn
6) tables_type: set to 3 for cubic lattice
7) taa1: nn hopping
8) tz: = taa1
9) V1a: disorder strength
10) liz: linear size of liz
11) The broadening eta is actually defined in the module_Global.for

3. The code initialize the liz lattice by reading the file lizList, which is created by a python script lizList.py.
In lizList.py, change the value of L(linear size of the cluster) and liz and run the python script

4. Running script:
In the current case, we set meas=5, run=1, so if we want to run a calculation with 400 disorder configuration,
we need to have 80 cores running in parallel, so that we can use the common
mpirun -np 80  ./AVEDISDCA_intel
A sample of running a system with L=4, liz=3 with CPA embedding is included in the folder sample.
In the folder, a submitting file sub.pbs on supermike is provided together with the output files.
The one we need is the dos.dat which contains the average dos(typical dos) for CPA(TMT) embedding case.
