#!/usr/bin/env python
# coding: utf-8


from collections import Counter
from itertools import combinations
from math import *
from random import *

import ase
import numpy as np
import pandas as pd
from pymatgen.core import Lattice, Structure


def gel_latt_coords(parameters):
    lattice = parameters[:6]
    coords = parameters[6:]
    coords = coords.reshape(int(coords.shape[0] / 3), 3)
    return lattice, coords


def LatticeObj(latticeparams, constrains):

    natoms = constrains["atoms"]
    vpa = constrains["vpa"]
    uh = (natoms * vpa[1]) ** (1 / 3)
    ul = (natoms * vpa[0]) ** (1 / 3)

    ub = np.array([uh for _ in range(3)] + [90 for _ in range(3)])
    lb = np.array([ul for _ in range(3)] + [90 for _ in range(3)])
    lattice = lb + (ub - lb) * latticeparams
    #    print(lattice)
    lattObj = Lattice.from_parameters(
        a=lattice[0],
        b=lattice[1],
        c=lattice[2],
        alpha=lattice[3],
        beta=lattice[4],
        gamma=lattice[5],
    )

    return lattObj


def createRandomData(constrains, trials=1000):

    natoms = constrains["atoms"]
    species = []
    for key in constrains["composition"].keys():
        nkey = int(
            (
                constrains["composition"][key]
                / sum(list(constrains["composition"].values()))
            )
            * natoms
        )
        for _ in range(nkey):
            species.append(key)

    count = 0
    while count < trials:
        shuffle(species)
        parameters = np.array([random() for _ in range(3 * natoms + 6)])
        structData = {"parameters": parameters, "species": species}
        if check_constrains(structData, constrains):
            break

    return structData


def check_constrains(structData, constrains):
    '''
    constrains = {
        "composition":{"Au":1},
        "atoms":13,
        "vpa":[16, 20],
        "r":r,
        }

    '''

    parameters = structData["parameters"].copy()
    species = structData["species"].copy()
    specieCount = dict(Counter(species))
    lattice, coords = gel_latt_coords(parameters)

    # print(species)

    natoms = constrains["atoms"]
    if natoms != coords.shape[0]:
        return False
    #    print(coords.flatten().tolist())

    composition = constrains["composition"]
    r = constrains["r"]

    latt = LatticeObj(lattice, constrains)

    M = np.array((latt.matrix))

    D = DistanceMatrix(coords, M)

    np.fill_diagonal(D, 1e300)
    DF = pd.DataFrame(D, columns=species, index=species)
    elements = r.columns.tolist()
    # print(DF)

    if len(list(composition.keys())) != len(list(specieCount.keys())):
        return False

    for key in composition.keys():
        if specieCount[key] % composition[key] != 0:
            return False

    for sp in elements:
        try:
            if (DF.loc[sp, sp].values < r.loc[sp, sp]).any():
                #               print("dist")
                return False
        except:
            if (DF.loc[sp, sp] < r.loc[sp, sp]).any():
                #                print("dist")
                return False

    for pair in combinations(elements, 2):
        try:
            if (
                DF.loc[pair[0], pair[1]].values < r.loc[pair[0], pair[1]]
            ).any():
                #                print("dist")
                return False
        except:
            if (DF.loc[pair[0], pair[1]] < r.loc[pair[0], pair[1]]).any():
                return False

    # print(structData)

    return True


def DistanceMatrix(frac_coordinates, M):

    a, b, c = (
        frac_coordinates[:, 0],
        frac_coordinates[:, 1],
        frac_coordinates[:, 2],
    )

    def getDist(mat):
        n, m = np.meshgrid(mat, mat)
        dist = m - n
        dist -= np.rint(dist)
        return dist

    da, db, dc = (
        getDist(a),
        getDist(b),
        getDist(c),
    )  # Fracrtional difference matrix

    # ---------cartesian differences------------

    DX = M[0][0] * da + M[1][0] * db + M[2][0] * dc
    DY = M[0][1] * da + M[1][1] * db + M[2][1] * dc
    DZ = M[0][2] * da + M[1][2] * db + M[2][2] * dc

    # -----------distance matrix--------------

    D = np.sqrt(np.square(DX) + np.square(DY) + np.square(DZ))

    return D


def SructureFrmParams(structData, constrains, pad=20):

    # print(structData ["parameters"])

    parameters = structData["parameters"].copy()
    species = structData["species"].copy()
    lattice, coords = gel_latt_coords(parameters)

    # print(lattice)

    lattice = LatticeObj(lattice, constrains)

    struct = Structure(lattice, species, coords, to_unit_cell=True)

    struct = add_padding(struct, pad=pad)

    return struct


def get_string(struct, energy):

    lattice = struct.lattice
    frac = (
        np.array([list(site.frac_coords) for site in struct.sites])
        .flatten()
        .tolist()
    )
    species = [site.specie.symbol for site in struct.sites]
    latt = [
        lattice.a,
        lattice.b,
        lattice.c,
        lattice.alpha,
        lattice.beta,
        lattice.gamma,
    ]

    dataString = (
        " ".join(map(str, latt + frac))
        + "|"
        + " ".join(species)
        + "|"
        + "{}".format(energy)
    )

    return dataString


def ParamsfromStruct(
    struct, constrains, energy=1e300, write_file="dumpfile.dat", pad=20
):

    dataString = get_string(struct, energy)

    with open(write_file, "a") as outfile:
        outfile.write("{}\n".format(dataString))

    struct = unpad(struct, pad=pad)

    lattice = struct.lattice
    frac = (
        np.array([list(site.frac_coords) for site in struct.sites])
        .flatten()
        .tolist()
    )
    species = [site.specie.symbol for site in struct.sites]
    latt = [
        lattice.a,
        lattice.b,
        lattice.c,
        lattice.alpha,
        lattice.beta,
        lattice.gamma,
    ]

    dataString = (
        " ".join(map(str, latt + frac))
        + "|"
        + " ".join(species)
        + "|"
        + "{}".format(energy)
    )

    latt = np.array(latt)

    natoms = constrains["atoms"]
    vpa = constrains["vpa"]
    uh = (natoms * vpa[1]) ** (1 / 3)
    ul = (natoms * vpa[0]) ** (1 / 3)

    latt[:3] = (latt[:3] - ul) / (uh - ul)
    latt[3:] = np.array([0, 0, 0])

    param = np.array(latt.tolist() + frac)

    structData = {"parameters": param, "species": species}

    return structData


def center(structure):

    center = np.average([s.frac_coords[2] for s in structure.sites])

    translation = (0.5 - center, 0.5 - center, 0.5 - center)

    structure.translate_sites(range(len(structure.sites)), translation)

    return structure


# --------------------------------


def add_padding(structure, pad=20):

    pos = [list(cord.frac_coords) for cord in structure.sites]

    species = [cord.specie.symbol for cord in structure.sites]
    latt = np.array(structure.lattice.matrix)

    multiplier_a = (pad + structure.lattice.a) / structure.lattice.a
    latt[0] = latt[0] * multiplier_a
    multiplier_b = (pad + structure.lattice.b) / structure.lattice.b
    latt[1] = latt[1] * multiplier_b
    multiplier_c = (pad + structure.lattice.c) / structure.lattice.c
    latt[2] = latt[2] * multiplier_c

    for i, site in enumerate(pos):
        pos[i][0] = site[0] / multiplier_a
        pos[i][1] = site[1] / multiplier_b
        pos[i][2] = site[2] / multiplier_c

    newstructure = Structure(latt, species, pos, to_unit_cell=True)
    return center(newstructure)


def unpad(structure, pad=20):

    pos = [list(cord.frac_coords) for cord in structure.sites]

    species = [cord.specie.symbol for cord in structure.sites]
    latt = np.array(structure.lattice.matrix)

    multiplier_a = (structure.lattice.a - pad) / structure.lattice.a
    latt[0] = latt[0] * multiplier_a
    multiplier_b = (structure.lattice.b - pad) / structure.lattice.b
    latt[1] = latt[1] * multiplier_b
    multiplier_c = (structure.lattice.c - pad) / structure.lattice.c
    latt[2] = latt[2] * multiplier_c

    for i, site in enumerate(pos):
        pos[i][0] = site[0] / multiplier_a
        pos[i][1] = site[1] / multiplier_b
        pos[i][2] = site[2] / multiplier_c

    newstructure = Structure(latt, species, pos, to_unit_cell=True)
    return center(newstructure)
