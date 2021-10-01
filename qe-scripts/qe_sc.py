from pymatgen.core.structure import IStructure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.pwscf import PWInput
import numpy as np
import argparse
import json
import sys
import os
import re

"""
Parsing arguments
"""

parser = argparse.ArgumentParser()
parser.add_argument('--tol', type=float, default='0.2', help='tolerance')
parser.add_argument('--path', type=str, default='.', help='Path to a folder with POSCAR')
parser.add_argument('--poscar', type=str, default='POSCAR', help='POSCAR file name')
parser.add_argument('--note', type=str, default='', help='note to add into script')
parser.add_argument('--o', type=bool, default=True, help='Overwrite all files')
parser.add_argument('--pressure', type=int, default=2000, help='Pressure in kBar')
args = parser.parse_args()

if args.path == '.':
    pwd = str(os.getcwd())
else:
    pwd = os.path.abspath(args.path)

"""
Reading structure, collecting space group and inequivalently positioned atoms
"""

structure = IStructure.from_file(os.path.join(pwd, args.poscar))
formula = structure.formula.replace(' ', '')
res = re.split('(\d+)', formula)
prefix = str()
short = str()
for c in res:
    if c.isalpha():
        prefix = prefix + c
    elif c.isnumeric():
        short = short + c
        if not c == '1':
            prefix = prefix + c
short = short + args.note

analyzer = SpacegroupAnalyzer(structure, symprec=args.tol)
sg = str(analyzer.get_space_group_number())

refined = analyzer.get_refined_structure()

a = round(refined.lattice.a, 10)
b = round(refined.lattice.b, 10)
c = round(refined.lattice.c, 10)
cosBC = round(np.cos(refined.lattice.alpha * np.pi / 180), 10)
cosAC = round(np.cos(refined.lattice.beta * np.pi / 180), 10)
cosAB = round(np.cos(refined.lattice.gamma * np.pi / 180), 10)

ref_analyzer = SpacegroupAnalyzer(refined, symprec=args.tol)
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

_kpoints = Kpoints.automatic_density(structure, kppa=50000)


def arr_str(s):
    return str(s).replace('[', '').replace(']', '').replace(',', '')


kpoints = arr_str(_kpoints.kpts[0])
shift = arr_str(_kpoints.kpts_shift)

"""
Writing input.opt
"""

input_opt = f"""&control
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
   press = {args.pressure},
   cell_dofree = 'xyz',
 /
ATOMIC_SPECIES
{species}
ATOMIC_POSITIONS (crystal_sg)
{positions}
K_POINTS {automatic}
 {kpoints}  {shift}""".format(automatic=automatic)

if args.o:
    with open(os.path.join(pwd, 'input.opt'), 'w') as f:
        f.write(input_opt)
        f.close()

"""
Writing script1.sh
"""

default_system = {
    "partition": "lenovo",
    "modules": "module load intel/mkl-11.2.3 mpi/impi-2018.2.199\nexport PATH=~/qe-6.6/bin/:s$PATH",
    "N_pw": 1,
    "n_pw": 20,
    "N_ph": 1,
    "n_ph": 40
}

if os.path.isfile(os.path.join(pwd, 'system.json')):
    with open(os.path.join(pwd, 'system.json'), 'r') as json_file:
        system = json.load(json_file)
        json_file.close()
else:
    system = default_system

script1 = f"""#!/bin/sh
#SBATCH -o qe.out1
#SBATCH -p {system['partition']}
#SBATCH -J 1L{short}
#SBATCH -N {system['N_pw']}
#SBATCH -n {system['n_pw']}

{system['modules']}

echo "OPT of {prefix} LAUNCHED at" $(date) | tee -a log.{prefix} 

mpirun $(which pw.x) -in $PWD/input.opt &> \$PWD/output.opt.${prefix} 

echo "OPT of {prefix} STOPPED at" $(date) | tee -a log.{prefix} 

""".format()

if args.o:
    with open(os.path.join(pwd, 'script1.sh'), 'w') as f:
        f.write(script1)
        f.close()