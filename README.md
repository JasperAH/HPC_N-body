#### Git repository contains:
- nbody implementation in Fortran90
- framework in python as provided by Kristian
- two attempts at tUPL implementations. The timesteps version seems usable for implementation in the python framework.

#### Todo:
- [x] compile and run Fortran90 implementation on a (relatively simple) dataset (nbody/input contains some)
    - run on the most simple (2_body.dat): `./bin/nbody input/2_body.dat 10 3 1 1e-6 data`
    - data gets output to nbody/data, containing final xyz positions and velocities and per-particle timestepped locations and velocities
    - using `gnuplot` and `plot.p`, the timestepped locations can be animated (might be nice for comparing to tUPL if using similar output format)
- [x] correct timestepping in timesteps.tUPL (currently wrong timestep gets updated at the end, maybe a separate pool like Kx,Kv,K_mask is necessary)
    - There are dedicated M, X and V shared spaces, X and V also require a time element to access, allowing for setting the correct X and V values in the 5th condition. N.B. this solution may be related to the issue described below since its cause is not certain for now.
- [x] implement tUPL in python framework: tuplcomplex.py seems suited for this as it has a similar structure
    - This is largely done, there is however an issue where not all orders of tupl execution finish. Some lead to deadlocks, not yet sure why for now. Might be to do with an attempt at forelem, since currently only whilelem is available.
    - [x] investigate & fix deadlock cause
- [x] add file reading such that both can use the same input data
- [x] add file writing to same format as Fortran90 format such that `gnuplot` script can be used
- [x] compare results of Fortran90 and tUPL implementation after running on the same dataset (THESE ARE VERY COMPARABLE /caps)
- [ ] rewrite and generalize if statements (both in tUPL spec and python) by adding only termination criteria for tupls instead of executing criteria
- [ ] add k-means functionality to script to cluster bodies
- [ ] add separate calculations for clusters and bodies
- [ ] integrate body and cluster calculations
- [ ] determine expected error (roughly, if necessary)
- [ ] tbd...
