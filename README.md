# Directions

You should be able to bring this repository into VSCode.  

In the terminal type "make" and it should make the build. 

Then inside "programs" directory you can run any program. Testing is just for testing certain functions. 

# Programs

## Gibbs Minimizer

Running the program gibbs_minimizer will create a file name "thermochemical_table.csv" that has 160001 lines. The first line is the header and shows the variables being written. 
This program uses a Gibbs Free-Energy minimization technique (using chemical potential and not the equilibrium constant method) to calculate thermodynamic variables to be used in 
my solve_real_gas program. It outputs density, internal energy, pressure, temperature, mixture gas constant, mizture c_v, mixture dp/drho, and mixture dp/de for an 10-species air model
(N2, O2, NO, N, O, Ar, Ar+, N+, O+, e-). It takes in a fixed density (rho) and internal energy (e) and uses NASA polynomials as well as iterative methods to minimize the free-energy. It
then takes the new species mass fractions and computes the above thermodynamic state variables from a newly solved equilibirum temperature that accounts for translational, rotational,
and vibrational energies as well as the enthalpy of formaton. There are two functions that you can swap between. One uses a constant density and varying internal energy to create a couple
of .csv files that can be used to plot the molar / mass fractions of each species for equilibrium temperatures ranging from ~400 K to 20,000 K. The second creates the thermochemical table.

## Perfect Gas Solver

You can run the program solve_perf_gas with up to two extra inputs. './solve_perf_gas' will run the hard coded grid, boundary conditions, and inlet conditions that you set in ../perfgassolver/main.cpp.
'./solve_perf_gas 400x200Ramp' will restart the file previously saved as '../plotfiles/PGI_400x200Ramp.dat'. You do not need to type in the 'PGI_' or '.dat', as this is already added. (PGI stands for 
Perfect Gas Inviscid, its not set up for viscous yet). Finally, './solve_perf_gas 400x200Ramp 50' will restart the 400x200 Ramp solution wwhile printing updates every 50 time steps. As stated, all of these
methods output a .dat file for plotting in Tecplot or Matlab. (Probably Paraview as well, but I haven't tried it.)

## Real Gas Solver

Running the program solve_real_gas will start solving the case set up in /realgassolver/main.cpp. It uses the thermochemical table that was created in the Gibbs minimizer program to solve an inviscid real-gas 
chemical equilibirum case. Outputs a .dat file just like the above. I am setting up a restart process like the above but it is not currently working.

