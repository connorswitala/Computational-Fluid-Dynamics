# Directions

You should be able to bring this repository into VSCode.  

In the terminal type "make" and it should make the build. 

Then inside "programs" directory you can run any program. Testing is just for testing certain functions. 

# Programs

## Gibbs Minimizer

Running the program gibbs_minimizer will create a file name "thermochemical_table.csv" that has 160001 lines. The first line is the header and shows the variabels being written. 
This program uses a Gibbs Free-Energy minimization technique to calculate thermodynamic variables to be used in my solve_real_gas program. It outputs density, internal energy, pressure,
temperature, mixture gas constant, mizture c_v, mixture dp/drho, and mixture dp/de for an 10-species air model (N2, O2, NO, N, O, Ar, Ar+, N+, O+, e-). It takes in a fixed density (rho) and
internal energy (e) and uses NASA polynomials as well as iterative methods to minimize the free-energy. It then takes the new species mass fractions and computes the above thermodynamic state variables 
from a newly solved equilibirum temperature that accounts for translational, rotational, and vibrational energies as well as the enthalpy of formaton. 

## Perfect Gas Solver

Running the program solve_perf_gas will bring up an option to use a preset case or set your own case. Entering "preset" will run the case that is set up in perfgassolver/main.cpp. Using "set" will run
you through a list of steps to set up the solver. It outputs a .dat file for plotting in Tecplot or Matlab. (Maybe Paraview as well, but I haven't tried it.)

## Real Gas Solver

Running the program solve_real_gas will start solving the case set up in /realgassolver/main.cpp. There is currently no set up configuration like the perfect gas solver has. It uses the thermochemical table 
that was created in the Gibbs minimizer program to solve an inviscid real-gas chemical equilibirum case. Outputs a .dat file just like the above. 

