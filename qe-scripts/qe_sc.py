from pymatgen.core.structure import IStructure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.pwscf import PWInput
import numpy as np
import argparse
import sys

"""
Parsing arguments
"""

parser = argparse.ArgumentParser()
parser.add_argument('--tol', type=float, default='0.2', help='tolerance')
parser.add_argument('--path', type=str, default='POSCAR', help='path to POSCAR file')


"""
Reading structure, collecting space group and inequivalently positioned atoms
"""

structure = IStructure.from_file('POSCAR_1-3-28_300')
prefix = structure.formula.replace(' ', '')
analyzer = SpacegroupAnalyzer(structure, symprec=0.2)
sg = str(analyzer.get_space_group_number())

refined = analyzer.get_refined_structure()
print(refined.lattice)

a = round(refined.lattice.a, 10)
b = round(refined.lattice.b, 10)
c = round(refined.lattice.c, 10)
cosBC = round(np.cos(refined.lattice.alpha*np.pi/180), 10)
cosAC = round(np.cos(refined.lattice.beta*np.pi/180), 10)
cosAB = round(np.cos(refined.lattice.gamma*np.pi/180), 10)

ref_analyzer = SpacegroupAnalyzer(refined, symprec=0.2)
symm_struc = ref_analyzer.get_symmetrized_structure()

nat = str(len(symm_struc.equivalent_sites))
ntyp = str(len(set(structure.species)))

species = list()

for specie in list(set(structure.species)):
      symbol = str(specie.symbol)
      mass = str(specie.atomic_mass).replace(' amu', '')
      upf = symbol + '.UPF'
      species.append(symbol + ' ' + mass + ' ' + upf)
species = "\n".join(species)

coords = list()

eq_sites = symm_struc.equivalent_sites
for site in eq_sites:
      specie = str(site[0].species)
      atom = ''.join([i for i in specie if not i.isdigit()])
      x = site[0].frac_coords[0]
      y = site[0].frac_coords[1]
      z = site[0].frac_coords[2]
      coords.append(atom + "\t" + "{:.10f}".format(x) + "\t" + "{:.10f}".format(y) + "\t" + "{:.10f}".format(z))
positions = "\n".join(coords)

automatic = '{automatic}'

"""
Writing input.opt
"""

input_opt = f"""
&control
    calculation='vc-relax'
    restart_mode='from_scratch',
    prefix={prefix},
    pseudo_dir = '.',
    outdir='.',
    wf_collect = .true.,
    nstep = 9999,
/
&system
        space_group={sg},
        A={a},
        B={b},
        C={c},
        cosAB={cosAB},
        cosBC={cosBC},
        cosAC={cosAC},
        nat={nat},
        ntyp={ntyp},
        ecutwfc=80.0,
        occupations='smearing',
        degauss=0.05,
 /
 &electrons
    conv_thr =  5.0d-7,
    mixing_beta = 0.6,
    electron_maxstep = 9999,
    diagonalization = 'cg',
 /
 &IONS
  ion_dynamics='bfgs',
 /
 &CELL
   cell_dynamics = 'bfgs',
   cell_factor = 2.5,
   press = 2000,
   cell_dofree = 'xyz',
 /
ATOMIC_SPECIES
{species}
ATOMIC_POSITIONS (crystal_sg)
{positions}
K_POINTS {automatic}
 24 24 24  0 0 0
""".format()

