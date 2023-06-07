import numpy as np
from ortools.linear_solver import pywraplp
from ortools.init import pywrapinit
import time

""" this grabs the known records of packing from a txt and puts them in an array. Not used yet. 
Find the txt files here : http://hydra.nat.uni-magdeburg.de/packing/csq/csq.html """
f = open("ratio.txt","r")
txt = f.readlines()
tabratio = np.zeros(10000)
for line in txt:
  n,fn = line.split()
  tabratio[int(n)-1]=float(fn)

f = open("radius.txt","r")
txt = f.readlines()
tabradius = np.zeros(10000)
for line in txt:
  n,fn = line.split()
  tabradius[int(n)-1]=float(fn)

""" Gets the number of grid points that would interfere with a given point. used to estimate M. """
""" Should be called with n an integer representing radius of circle if grid was of size 1 """
""" Call it with n = R * delta^-1 """
#source : https://oeis.org/A000328
def Pointsincircle(n):
    return (sum([int((n**2 - y**2)**0.5) for y in range(1, n)]) * 4 + 4*n + 1)

""" Testing. Ignore """
"""testdelta1 = 10 #1/delta
testdelta = 1/testdelta1
testR = 0.2 

M1 = Pointsincircle(int(2/testdelta))+1
M2 = int(np.pi*((4/(testdelta**2))+(np.sqrt(2)/testdelta))) 
print(M1,M2)

print(4 * (testR*testdelta1)**2)"""
""" """ """ """ """ """
def bestsolutionforR(R):
    i=0
    while tabradius[i] > R:
        i+=1
    return i


def interference(i1,j1,i2,j2,R,dens):
    """ circles of radius R, in a grid of density d (d = 1/delta), with centers (i1,j1) and (i2,j2), will overlap when this condition is met"""
    return (i1 - i2)**2 + (j1 - j2)**2 < (4 * (R*dens)**2)



def pack(R,density):
    """ creates a solver """
    # IMPORTANT # If you don't have gurobi set up, replace GUROBI_MIP with SCIP and comment the "LoadGurobiSharedLibrary" line at the very end
    solver = pywraplp.Solver.CreateSolver('GUROBI_MIP')
    if not solver:
        return
    
    N = int((1-2*R)*density) # number of points per row/column

    """ create an array of variables """
    x = {}
    for i in range(N):
        for j in range(N):
            x[(i, j)] = solver.IntVar(0, 1, 'x_%i_%i' % (i, j))

    

    """ create constraints """
    M = Pointsincircle(int(R*density)+1) * 4
    print(M)
    for i in range(N):
        #print(i)
        for j in range(N):
            constraint = solver.RowConstraint(0, M, '')
            constraint.SetCoefficient(x[(i,j)], M)
            for i2 in range(N):
                for j2 in range(N):
                    if(interference(i,j,i2,j2,R,density) and (i != i2 or j != j2)): #only points that are in the neighborhood of (i,j) appear in the constraint
                        constraint.SetCoefficient(x[(i2,j2)], 1)

    print('Number of variables =', solver.NumVariables())
    print('Number of constraints =', solver.NumConstraints())


    starttime = time.time()
    solver.Maximize(solver.Sum([solver.Sum([x[(i,j)] for j in range(N)]) for i in range(N)]))
    status = solver.Solve()
    timetaken = time.time() - starttime

    if status == pywraplp.Solver.OPTIMAL:
        print('Solution:')
        print('Objective value =', solver.Objective().Value())
        for v in solver.variables():
            if v.solution_value() == 1:
                print(v.name())
            elif v.solution_value() != 0:
                print("ERROR", v.name, v.solution_value())   

        print('Objective value =', solver.Objective().Value()) 
        print('Best known number of circles of this radius : %d' % bestsolutionforR(R))
    else:
        print('The problem does not have an optimal solution.')

    
    print('\nAdvanced usage:')
    print('Problem solved in %f seconds' % timetaken)
    print('Problem solved in %d iterations' % solver.iterations())
    print('Problem solved in %d branch-and-bound nodes' % solver.nodes()) 
    




def main():
    """ "density" is to be understood as 1/delta, with delta being the gap between neighbouring points of the grid. """
    """ "R" is the radius of the circles. """
    Radius = 0.1
    pack(R=Radius,density=25)
    

if __name__ == '__main__':

    pywrapinit.CppBridge.InitLogging('intoarray.py')
    cpp_flags = pywrapinit.CppFlags()
    cpp_flags.logtostderr = True
    cpp_flags.log_prefix = False
    pywrapinit.CppBridge.SetFlags(cpp_flags)
    pywrapinit.CppBridge.LoadGurobiSharedLibrary('C:\gurobi952\win64\lib\gurobi95.lib')
    main()