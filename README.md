# PyAS15 (PPPL23)


## Description

The goal of this project is to explain some lower order calculations based on the IAS15 algorithm from [Rein and Spiegel 2015](https://arxiv.org/abs/1409.4779) and implemented in the [REBOUND code](https://rebound.readthedocs.io/en/latest/integrators/). The aim is for readability and understanding to explain this algorithm. The complexity is as follows:

The tests call functions from these files.
- ias15FirstOrder1D.py - equations of the form dy/dx = f(y)
- ias15FirstOrder2D.pyp - equations of the form dy/dx = f(y,x)
- ias15PDE.py - Equations of hamiltonian form dy/dt = f(x,y) and dx/dt = g(x,y)
- rk4.py - a simple rk4 integrator to compare against
- euler.py - a simple euler integrator to compare against.
Regression tests compare against an exponential growth problem and tracing circular orbits in phase space.


## Installation

Using conda, one can run these files with just numpy and matplotlib.pyplot.


## Contact

- [GitHub](https://github.com/maxwell-rosen)
- [Email](mailto:mhrosen@pppl.gov)
