import numpy as np
from matplotlib import pyplot as plt  
from PIL import Image, ImageDraw
from ortools.linear_solver import pywraplp
from ortools.init import pywrapinit
import time

""" this grabs the known records of packing from a txt and puts them in an array.
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
    """ grabs the best known packing of circles of radius R, using the table generated at the beginning"""
    i=0
    while (tabradius[i] >= R or tabradius[i] ==0):
        i+=1
    return i


def interference(i1,j1,i2,j2,R,dens):
    """ circles of radius R, in a grid of density d (d = 1/delta), with centers (i1,j1) and (i2,j2), will overlap when this condition is met"""
    return (i1 - i2)**2 + (j1 - j2)**2 < (4 * (R*dens)**2)



def pack(R:float,density:int,text:bool,s:int):
    """ Packs as many circles of radius R with grid of given density. set text to true for information ; then s decides the size of the image. s=0 generates no image """
    
    
    """ creates a solver """
    # IMPORTANT # If you don't have gurobi set up, replace GUROBI_MIP with SCIP and comment the "LoadGurobiSharedLibrary" line at the very end
    """ Gurobi is MUCH faster, but requires a license. """
    solver = pywraplp.Solver.CreateSolver('GUROBI_MIP')
    if not solver:
        return
    

    N = int((1-2*R)*density) + 1 # number of points per row/column

    """ create an array of variables """
    x = {}
    for i in range(N):
        for j in range(N):
            x[(i, j)] = solver.IntVar(0, 1, 'x_%i_%i' % (i, j))

    

    """ create constraints """
    M = Pointsincircle(int(R*density)+1)
    if(text):
        print("M : ", M)
        """ the constraints are of the form : M(1-x_ij) >= sum_(i',j' s.t. x_ij and x_i'j' interfere) x_i'j' """
        print("number of points per row : ", N)
    for i in range(N):
        for j in range(N):
            constraint = solver.RowConstraint(0, M, '')
            constraint.SetCoefficient(x[(i,j)], M) 
            for i2 in range(N):
                for j2 in range(N):
                    if(interference(i,j,i2,j2,R,density) and (i != i2 or j != j2)): #only points that are in the neighborhood of (i,j) appear in the constraint
                        constraint.SetCoefficient(x[(i2,j2)], 1)

    if(text):
        print('Number of variables =', solver.NumVariables())
        print('Number of constraints =', solver.NumConstraints())

    print("calling solver...")
    starttime = time.time()
    solver.Maximize(solver.Sum([solver.Sum([x[(i,j)] for j in range(N)]) for i in range(N)]))
    status = solver.Solve()
    timetaken = time.time() - starttime

    if status == pywraplp.Solver.OPTIMAL:
        bestknown = bestsolutionforR(R)
        if(text):
            print('Solution:')
            print('Objective value =', solver.Objective().Value())
            """ 
            for v in solver.variables():
                if v.solution_value() == 1:
                    print(v.name())
                elif v.solution_value() != 0:
                    print("ERROR", v.name, v.solution_value())   
            """
            print('Objective value =', solver.Objective().Value()) 
            
            print('Best known number of circles of this radius : %d' % bestknown)
            
            if(s > 0):
                """ Draws the actual image. s is sizeof the image in pixels. """
                img = Image.new('RGB', (s,s))
                d = ImageDraw.Draw(img)
                delta = 1 / density
                for i in range(N):
                    for j in range(N):
                        if x[(i, j)].solution_value() == 1:

                            d.point([s*((i*delta)+R), s*((j*delta)+R)],'Red')
                            d.ellipse([s*((i*delta)),s*((j*delta)),s*((i*delta)+2*R),s*((j*delta)+2*R)])
                        else:

                            d.point([s*((i*delta)+R), s*((j*delta)+R)],'Blue')
                img.show('Radius : %f ; Precision : %d' % (R,density) )
            else:
                """ prints a simple grid """
                str = ''
                for i in range(N):
                    for j in range(N):
                        if x[(i, j)].solution_value() == 1:
                            str += 'O'
                        else:
                            str += '·' 
                    str += '\n'
                print(str)
    else:
        print('The problem does not have an optimal solution.')

    
    
    print('Problem solved in %f seconds' % timetaken)
    print('Problem solved in %d iterations' % solver.iterations())
    print('Problem solved in %d branch-and-bound nodes' % solver.nodes()) 
    ratio = solver.Objective().Value()/bestknown
    return(ratio, timetaken)




def main():
    """ "density" is to be understood as 1/delta, with delta being the gap between neighbouring points of the grid. """
    """ "R" is the radius of the circles. """
    
    """ This generates a graph of ratios and time taken for the best radiuses for 1 to 100 circles """
    """ expratio = np.zeros(60)
    exptime = np.zeros(60)
    for k in range(1,60):
        Radius = 0.009
        print(k)
        expratio[k], exptime[k] = pack(R=Radius,density=k,text=False,s=0)
    plt.subplot(2, 1, 1)
    plt.ylim(0,1)
    plt.plot(expratio)
    plt.title('Ratio')  
    plt.subplot(2, 1, 2) 
    plt.plot(exptime)
    plt.title('Solver time')
    plt.show()    """
    pack(0.04,50,True,4000)
    

    
        
    

if __name__ == '__main__':

    pywrapinit.CppBridge.InitLogging('intoarray.py')
    cpp_flags = pywrapinit.CppFlags()
    cpp_flags.logtostderr = True
    cpp_flags.log_prefix = False
    pywrapinit.CppBridge.SetFlags(cpp_flags)
    #put the path of gurobi library here, or comment this next line if you don't have one.
    pywrapinit.CppBridge.LoadGurobiSharedLibrary('C:\gurobi952\win64\lib\gurobi95.lib')
    main()