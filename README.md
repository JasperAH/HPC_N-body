#### Git repository contains:
- nbody implementation in Fortran90
- framework in python as provided by Kristian
- two attempts at tUPL implementations. The timesteps version seems usable for implementation in the python framework.

#### Todo:
- compile and run Fortran90 implementation on a (relatively simple) dataset (nbody/input contains some)
- correct timestepping in timesteps.tUPL (currently wrong timestep gets updated at the end, maybe a separate pool like Kx,Kv,K_mask is necessary)
- implement tUPL in python framework: tuplcomplex.py seems suited for this as it has a similar structure
- add file reading/writing such that both can use the same input data
- compare results of Fortran90 and tUPL implementation after running on the same dataset
- tbd...
