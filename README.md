# Ising2025
Multi-threading simulation of the 2d Ising model with an external magnetic field.

The C++ code computes the numerical simulations of the 2d Ising model. The Python code generates a movie of the system dynamics.</br></br>
<b>Exportations:</b> spin-state snapshots and time-evolution of the total magnetization.</br>
<b>Compile:</b> g++ Ising_omp.cpp -fopenmp -lgsl -lgslcblas -lm -O3 -s -o Ising_omp.out.</br>
<b>Run:</b> ./Ising_omp.out -parameter=value.</br>
<b>Generate the movie:</b> python figure_Ising_dynamics.py -parameter=value.</br>
<b>List of parameters:</b> beta, h, LX, LY, tmax, init, ran, threads (details as comments in the code).
