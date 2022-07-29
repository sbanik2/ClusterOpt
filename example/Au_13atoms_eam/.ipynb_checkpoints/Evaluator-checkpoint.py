#!/usr/bin/env python
# coding: utf-8

# In[1]:


import math

from lammps import lammps
from pymatgen.io.lammps.data import LammpsData

from Structfunc import ParamsfromStruct, SructureFrmParams, check_constrains

# In[ ]:

cmds = ["-screen", "log.screen"]
lmp = lammps(cmdargs=cmds)
# lmp  = lammps()


class Custom_minimize(object):
    def __init__(
        self,
        fun,
        x0,
        args,
        jac=None,
        hess=None,
        hessp=None,
        bounds=None,
        constraints=None,
        callback=None,
        **options
    ):

        self.parameters = x0
        self.constrains = args[0]
        self.lammps_args = args[1]
        self.pair_style = self.lammps_args["pair_style"]
        self.pair_coeff = self.lammps_args["pair_coeff"]
        self.pad = self.lammps_args["pad"]
        
        self.success = True
        
        try:
            self.x, self.fun = self.minimize()
        except:
            self.success = False
            
            
    def minimize(self):

        structData = {
            "parameters": self.parameters,
            "species": [
                list(self.constrains["composition"].keys())[0] for _ in range(self.constrains["atoms"])
            ],
        }
        
        #print(structData ["species"])
        #print(structData ["parameters"].copy())
        if not check_constrains(structData, self.constrains):
            print("constraind failed")
            return  self.parameters,1e300

        struct = SructureFrmParams(
            structData, self.constrains, pad=self.pad
        )
        LammpsData.from_structure(struct, atom_style="atomic").write_file(
            "in.data"
        )

        lmp.command("clear")
        lmp.command("dimension 3")
        lmp.command("box tilt large")
        lmp.command("units metal")
        lmp.command("atom_style atomic")
        lmp.command("neighbor 2.0 bin")
        lmp.command("atom_modify map array sort 0 0")
        lmp.command("boundary f f f")
        lmp.command("read_data in.data")
        lmp.command("{}".format(self.pair_style))
        lmp.command("{}".format(self.pair_coeff))
        lmp.command("thermo 1000")
        lmp.command("thermo_style custom step etotal atoms vol")
        lmp.command("thermo_modify format float %5.14g")
        lmp.command("variable potential equal pe/atoms")
        lmp.command("neigh_modify one 5000 delay 0 every 1 check yes")
        lmp.command("run 0 pre no")
        tmp_eng = lmp.extract_variable('potential', None, 0)
        if math.isnan(tmp_eng) or math.isinf(tmp_eng):
            lmp.command("write_data min.geo")
        else:
            lmp.command("minimize 1.0e-8 1.0e-8 10000 10000")
            lmp.command("write_data min.geo")

        lmp.command("run 0 pre no")

        energy = lmp.extract_variable('potential', None, 0)

        if math.isinf(float(energy)):
            energy = 1e300
        if math.isnan(float(energy)):
            energy = 1e300

        minStruct = LammpsData.from_file(
            "min.geo", atom_style="atomic"
        ).structure
        outData = ParamsfromStruct(
            minStruct,
            self.constrains,
            energy=energy,
            write_file="dumpfile.dat",
            pad=self.pad,
        )
        return outData ["parameters"],energy 


class LammpsEvaluator(object):
    '''
    Performs energy evaluations on a offspring
    '''

    def __init__(
        self,
        constrains,
        lammps_args,
    ):

        self.constrains = constrains
        self.lammps_args = lammps_args
        self.pair_style = self.lammps_args["pair_style"]
        self.pair_coeff = self.lammps_args["pair_coeff"]
        self.pad = self.lammps_args["pad"]

    def snapshot(self, parameters):

        structData = {
            "parameters": parameters,
            "species":  [
                list(self.constrains["composition"].keys())[0] for _ in range(self.constrains["atoms"])
            ],
        }

        if not check_constrains(structData, self.constrains):
            return 1e300

        struct = SructureFrmParams(structData, self.constrains, pad=self.pad)
        LammpsData.from_structure(struct, atom_style="atomic").write_file(
            "in.data"
        )

        lmp.command("clear")
        lmp.command("dimension 3")
        lmp.command("box tilt large")
        lmp.command("units metal")
        lmp.command("atom_style atomic")
        lmp.command("neighbor 2.0 bin")
        lmp.command("atom_modify map array sort 0 0")
        lmp.command("boundary f f f")
        lmp.command("read_data in.data")
        lmp.command("{}".format(self.pair_style))
        lmp.command("{}".format(self.pair_coeff))
        lmp.command("thermo 1000")
        lmp.command("thermo_style custom step etotal atoms vol")
        lmp.command("thermo_modify format float %5.14g")
        lmp.command("variable potential equal pe/atoms")
        lmp.command("neigh_modify one 5000 delay 0 every 1 check yes")
        lmp.command("run 0 pre no")
        energy = lmp.extract_variable('potential', None, 0)
        if math.isinf(float(energy)):
            energy = 1e300
        if math.isnan(float(energy)):
            energy = 1e300

        return energy
