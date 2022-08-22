# ClusterOpt

<p align="justify"> Implementation of Basin Hopping Global Optimization method for optimization of Atomic nanoclusters. </p>

The following paper describes the details of the CGCNN framework:

## Table of Contents
- [Introduction](#Introduction)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Running the code](#Running-the-code)
- [Citation](#data-availability)
- [License](#license)

## Introduction
 This package implements the Basin Hopping global optimization method in scipy (https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.basinhopping.html) along with LAMMPS (https://www.lammps.org/) MD (Molecular dynamic simulation) package to find the global minima of atomic system (Single or Multicomponent)

<p align="center"> <a href="url"><img src="https://github.com/sbanik2/CEGAN/blob/main/Figs/Workflow.png" align="center" height="400" width="600" ></a> </p>



## Prerequisites
This package requires:
- [scipy](https://scipy.org/)
- [LAMMPS](https://www.lammps.org/)
- [pymatgen](https://pymatgen.org/)
- [pandas](https://pandas.pydata.org/)
- [numpy](https://numpy.org/)


## Installation

### Manual Installation
Install the anaconda package (https://docs.anaconda.com/anaconda/install/). Then, 
```
conda env create --name ClusterOpt -f environment.yml
conda activate ClusterOpt
git clone https://github.com/sbanik2/ClusterOpt.git
python setup.py install
```

### Installation with pypi

```
pip install ClusterOpt

```
<p align="justify"> ***The package requires lammps binding to run. First lammps package needs to be downloaded from  https://www.lammps.org/download.html and compiled. The instructions on python integration can be found here https://docs.lammps.org/Python_install.html.</p>


### Running the code
<p align="justify"> An example of running run directory provided in the example section. First all the parameters crystal and the lammps pair_style and pair_coeff should be set. The composition is given for e.g., a Au<2>Al<3> as "composition":{"Au":2,"Al:3"}, the minimum interatomic distances as a pandas dataframe with rows and columns belonging to each species in the same order they are mentioned in the composition. E.g.,</p>

``` python
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

```

Once the parameters are set, the LammpsEvaluator can be initialized and the basinhopping optimizer can be set for the run

 ``` python

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
```
details of individual hyperparameters for the basinhopping optimizer can be found here https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.basinhopping.html.



### Citation
```
@article{banik2022cegan,
  title={CEGAN: Crystal Edge Graph Attention Network for multiscale classification of materials environment},
  author={Banik, Suvo and Dhabal, Debdas and Chan, Henry and Manna, Sukriti and Cherukara, Mathew and Molinero, Valeria and Sankaranarayanan, Subramanian KRS},
  journal={arXiv preprint arXiv:2207.10168},
  year={2022}
}
```
### License
CEGAN is licensed under the MIT License
![image](https://user-images.githubusercontent.com/66140668/185830167-753fbfcd-76fa-4c55-8140-41216a2a713d.png)

