import pandas as pd
import numpy as np
from random import *
from scipy.optimize import minimize
from Evaluator import Custom_minimize, LammpsEvaluator
from Structfunc import createRandomData
from scipy.optimize import basinhopping


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
             niter=10000000,
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
#             target_accept_rate=0.5, 
#            stepwise_factor=0.9
            )