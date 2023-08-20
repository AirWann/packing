# An implementation of grid-based packing discretization

This projects uses GLOP, Google's Linear Optimization Package, with the Gurobi solver, to implement the techniques explained in the internship report. 
If you do not have a Gurobi licence, some lines have to be changed for the code to work. There are comments clearly stating which.

## Reading images

Plots shown in the images folder were generated using this code.
"Ratio" means the ratio of approximation of the solution, as in, the number of circles packed divided by the best known number of circles packed. these best known solutions were grabbed from http://hydra.nat.uni-magdeburg.de/packing/csq/csq.html. The time, in seconds, is only for _solving_ the LP and not generating it.

