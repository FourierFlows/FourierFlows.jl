# FourierFlows

This software provides a framework for solving doubly-periodic problems 
using Fourier-based pseudospectral methods. The code is written in Julia.

To provide a convenient a framework as possible for rapidly developing 
fast, optimized solvers for two-dimensional problems, the code is divided
along conceptual lines into problem-agnostic and problem-specific components. 
The problem-agnostic parts of the code are contained in ``/src/``. The problem-
specific parts of the code are single files, one for each problem, in the directory
``/src/physics``.

# Authors

Gregory L. Wagner and Navid C. Constantinou
