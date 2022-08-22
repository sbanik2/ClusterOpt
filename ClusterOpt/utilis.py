import os
import sys
import numpy as np
from pymatgen.io.vasp.inputs import Poscar
from ClusterOpt.Structfunc import gel_latt_coords 
from pymatgen.core import Structure,Lattice
import math




def CreateStructure(filename,path,nStructure=1):

    if not os.path.exists(path):
        os.mkdir(path)

    structurelist = []

    
    with open(filename,"r") as infile:
        for i,line in enumerate(infile):
            col = line.split("|")

            parms = np.array([float(val) for val in col[0].split()])
            species = col[1].split()
            energy = float(col[2])

            if math.isnan(energy):
                continue


            structurelist.append([energy,parms,species])

    structurelist.sort(key=lambda x:x[0])


    for i,val in enumerate(structurelist[:nStructure]):

        parameters = val[1]
        species = val[2]
        lattice,coords = gel_latt_coords(parameters)


        lattice = Lattice.from_parameters(a=lattice[0],b=lattice[1],
                                          c=lattice[2],alpha=lattice[3],
                                          beta=lattice[4],gamma=lattice[5])



        struct = Structure(lattice, species ,coords,to_unit_cell=True)


        with open("energy.dat","a") as outfile:
            outfile.write("{} {}\n".format(i,val[0]))

        Poscar(struct).write_file("{}/{}.POSCAR".format(path,i))

        print("structure {} created in {}".format(i,path))
        

        
def Status(outfile):
    
    with open(outfile, "r") as infile:
        minval = 1e7
        coords = ""
        cnt = 0
        minloc = 0
        for i, line in enumerate(infile):
            cnt += 1
            col = line.split("|")
    #        col = line.split("|")
            eng = float(col[-1])
            if minval > eng:
                minval = eng
                minloc = i+1
                coords = col[1]
    
    print("Lowest Value: %s"%(minval))
    print("Number of Total Evaluations: %s"%(cnt))
    print("Number of Evaluations till Minima was found: %s"%(minloc))


