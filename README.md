# STG
A project for testing varoius methods of synthetic turbulence generation (STG methods in further). 
Such methods are employed in eddy-resolved methods of Navier-Stoks equations solution like LES, DNS etc. 
They are used in such methods mainly for two purposes: 
1. Imposing of turbulence inlet boundary conditions 
2. For imposing velocity pulsation at the inlet of LES zone in RANS-LES intefaces in Embedded LES methods. 

Five different STG methods were realised in generally:
1. Method of "white noise" generation proposed by Lund in [1]
2. Spectral method of Smirnov proposed in [2]
3. Spectral method of Davidson proposed in [3]
4. Synthetic eddy method (SEM) of Jarrin prposed in [4]
5. Spectral method based on Davidson's method but with explicit dependency of velocity from time and with 
capability for user to specify velocity spectrum.

### Running tests
Methods were realised in C. C files are contained in *STG_core* directory. 
To build shared library from these files you need to run *build.sh* from project directory. 
After completing of *build.sh* working *build* folder in project directory will appear. 
Library files will be placed in it.

Tests were wrote in Python with utilization of ctypes library for calling function from shared library with STG methods. 
Python project with test locates in *STG_test_tool* directory. 
To run tests user need to run command ```pytest tests.py``` from STG_test_tool directory. 
To get more information about invocation of tests created 
by means of pytest library you can look [pytest documentation](http://doc.pytest.org/en/latest/usage.html). 
After tests completing *plots* folder will appear in *STG_test_tool* directory. 
Plots of different statistical parameters of generated velocity field will saved in this folder.

### Configuring tests
For testing of STG methods generation of velocity on uniform 2D mesh on a given time interval is conducted. 
Then different statistical parameters like second moments in point or correlation (space and time) 
functions are computed.
Velocity field is homogeneous. 
You can specify parameters of the field (Reynolds stresses tensor, time and length scales) 
in *config.py* file. Also in this file settings of time intervals and mesh are specified. 
Settings of plots (line width, axes limits etc.) are specified in *config_plots.py* file.
Methods specific parameters like number for modes in spectral methods or number of eddies in SEM are specified in tests 
declarations in *tests.py* file.

### Requirements
1. Anaconda 3 with Python of version 3.7.0 and higher
2. GCC compiler. (Version 7.4.0 was used for developing)
## References
**[1]** T. Lund, X. Wu, and D. Squires Generation of turbulent inflow data for spatially-developing boundary 
        layer simulations // Journal of Computational Physics, 1998, N. 140, P. 233–258.
 
**[2]** A. Smirnov, S. Shi, and I.B. Celik Random flow generation technique for large-eddy simulations and 
        particle-dynamics modeling // Journal of Fluids Engineering, 2001, N. 123, P. 359–371.

**[3]** L. Davidson Using isotropic synthetic fluctuations as inlet boundary conditions for unsteady 
        simulations // Advances and Applications in Fluid Mechanics, 2007, N. 1, P. 1–35.

**[4]** N. Jarrin, S. Benhamadouche, D. Laurence, and R. Prosser,A synthetic-eddy-method for generating inflow 
        conditions for large-eddy simulations // International 
        Journal of Heat and Fluid Flow, 2006, V. 27, pp. 585–593.


