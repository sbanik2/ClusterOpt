#ClusterOpt

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
<p align="justify">  This package implements the Basin Hopping global optimization method in scipy(https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.basinhopping.html.) along with LAMMPS (https://www.lammps.org/#gsc.tab=0) MD (Molecular dynamic simulation) package to find the global minima of atomic system (Single or Multicomponent) </p>

<p align="center"> <a href="url"><img src="https://github.com/sbanik2/CEGAN/blob/main/Figs/Workflow.png" align="center" height="400" width="600" ></a> </p>



## Prerequisites
This package requires:
- [scipy](https://scipy.org/)
- [LAMMPS](https://www.lammps.org/)
- [pymatgen](https://pymatgen.org/)
- [pandas](https://pandas.pydata.org/)
- [numpy](https://numpy.org/)


## Installation

#Manual Installation
Install the anaconda package (https://docs.anaconda.com/anaconda/install/). Then, 
```
conda env create --name ClusterOpt -f environment.yml
conda activate ClusterOpt
```

For installation with pypi

```
pip install ClusterOpt

```
***The package requires lammps binding to run. First lammps package needs to be downloaded from and compiled and the instructions https://www.lammps.org/download.html on python integration can be found here https://docs.lammps.org/Python_install.html.


### Running the code
An example of running run directory provided in the example section. First all the parameters crystal and the lammps pair_style and pair_coeff should be set

.. code:: python
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

