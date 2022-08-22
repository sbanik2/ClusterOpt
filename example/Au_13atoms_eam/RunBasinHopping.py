import numpy as np
import pandas as pd
from ClusterOpt.Evaluator import Custom_minimize, LammpsEvaluator
from ClusterOpt.Structfunc import createRandomData
from ClusterOpt.utilis import CreateStructure,Status
from scipy.optimize import basinhopping
from scipy.optimize import minimize


# minimum distance criteria between the atoms

r = pd.DataFrame(np.array([[1.3]]),columns=["Au"],index= ["Au"])

constrains = {     
        "composition":{"Au":1},
        "atoms":13,
        "vpa":[16, 20],
        "r":r,       
        }


lammps_args = {
        "pair_style" : "pair_style eam",
        "pair_coeff" : "pair_coeff * * Au.eam",
        "pad":20
            }

args = (constrains,lammps_args)


structData = createRandomData(constrains,trials = 1000)
x0 = structData["parameters"]


fun = LammpsEvaluator(constrains,lammps_args).snapshot
res = minimize(fun, x0, args=args, method=Custom_minimize)

print(res.x,res.fun)

minimizer_kwargs = {"method":Custom_minimize,"args":args}

def print_fun(x, f, accepted):
    print("at minimum %.2e accepted %d" % (f, int(accepted)))

basinhopping(fun,
             x0, 
             niter=1000,
             T=1.0, 
             stepsize=0.5, 
             minimizer_kwargs=minimizer_kwargs, 
             take_step=None, 
             accept_test=None, 
             callback=print_fun, 
             interval=50, 
             disp=False, 
             niter_success=None, 
             seed=None, 
            )

Status("dumpfile.dat")


CreateStructure("dumpfile.dat","./",nStructure=1)



