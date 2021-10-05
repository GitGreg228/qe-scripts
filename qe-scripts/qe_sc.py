from pymatgen.core.structure import IStructure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.inputs import Kpoints
import numpy as np
import argparse
import shutil
import json
import os
import re

"""
Parsing arguments
"""


def boolean_string(s):
    if s not in {'False', 'True'}:
        raise ValueError('Not a valid boolean string')
    return s == 'True'


parser = argparse.ArgumentParser()
parser.add_argument('--tol', type=float, default='0.2', help='tolerance')
parser.add_argument('--path', type=str, default='.', help='Path to a folder with POSCAR')
parser.add_argument('--poscar', type=str, default='POSCAR', help='POSCAR file name')
parser.add_argument('--note', type=str, default='', help='note to add into script')
parser.add_argument('--o', type=bool, default=True, help='Overwrite all files')
parser.add_argument('--press', nargs='+', default=2000, help='Pressure(s) in kBar')
parser.add_argument('--kppa', type=int, default=50000, help='K-points per unit volume')
parser.add_argument('--primitive', default=True, type=boolean_string, help='if print primitive structure')
parser.add_argument('--dyn', type=str, default='ibrav', help='cell dynamics')
args = parser.parse_args()

if args.path == '.':
    pwd = str(os.getcwd())
else:
    pwd = os.path.abspath(args.path)


def arr_str(s):
    return str(s).replace('[', '').replace(']', '').replace(',', '').replace('\'', '')


def overwrite(path, name, content, o):
    tmp_path = os.path.join(path, name)
    if not os.path.isfile(tmp_path):
        with open(tmp_path, 'w') as f:
            f.write(content)
            f.close()
    elif o:
        with open(tmp_path, 'w') as f:
            f.write(content)
            f.close()


def formulas(structure):
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
    return prefix, short


def create_input_opt(tol, path, structure, note, o, pressure, kppa, dyn):
    """
    Reading structure, collecting space group and inequivalently positioned atoms
    """

    prefix = formulas(structure)[0]
    short = formulas(structure)[1]
    short = short + note
    analyzer = SpacegroupAnalyzer(structure, symprec=tol)
    sg = str(analyzer.get_space_group_number())

    refined = analyzer.get_refined_structure()

    a = round(refined.lattice.a, 10)
    b = round(refined.lattice.b, 10)
    c = round(refined.lattice.c, 10)
    cosBC = round(np.cos(refined.lattice.alpha * np.pi / 180), 10)
    cosAC = round(np.cos(refined.lattice.beta * np.pi / 180), 10)
    cosAB = round(np.cos(refined.lattice.gamma * np.pi / 180), 10)

    ref_analyzer = SpacegroupAnalyzer(refined, symprec=tol)
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

    _kpoints = Kpoints.automatic_density(structure, kppa=kppa)

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
        degauss=0.015,
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
   press = {pressure},
   cell_dofree = '{dyn}',
 /
ATOMIC_SPECIES
{species}
ATOMIC_POSITIONS (crystal_sg)
{positions}
K_POINTS {automatic}
 {kpoints}  {shift}""".format(automatic=automatic)

    overwrite(path, 'input.opt', input_opt, o)

    """
    Writing script1.sh
    """

    default_system = {
        "partition": "lenovo",
        "modules": "module load intel/mkl-11.2.3 mpi/impi-2018.2.199\nexport PATH=~/qe-6.6/bin/:s$PATH",
        "N_pw": 1,
        "n_pw": 20,
        "N_ph": 1,
        "n_ph": 40,
        "pp_path": "../PP"
    }

    if os.path.isfile(os.path.join(path, 'system.json')):
        with open(os.path.join(path, 'system.json'), 'r') as json_file:
            system = json.load(json_file)
            json_file.close()
    else:
        system = default_system

    if os.path.isdir(system['pp_path']):
        for fname in os.listdir(system['pp_path']):
            if '.UPF' in fname:
                tmp_src = os.path.join(system['pp_path'], fname)
                if os.path.isfile(tmp_src):
                    shutil.copy2(tmp_src, path)

    script1 = f"""#!/bin/sh
#SBATCH -o qe.out1 -e qe.err1
#SBATCH -p {system['partition']}
#SBATCH -J r{short}
#SBATCH -N {system['N_pw']}
#SBATCH -n {system['n_pw']}

{system['modules']}

echo "OPT of {prefix} LAUNCHED at" $(date) | tee -a log.{prefix} 

$(which mpirun) $(which pw.x) -in $PWD/input.opt &> $PWD/output.opt.{prefix} 

echo "OPT of {prefix} STOPPED at" $(date) | tee -a log.{prefix} 

""".format()

    overwrite(path, 'script1.sh', script1, o)


print('Reading structure from ' + os.path.join(pwd, args.poscar))
structure = IStructure.from_file(os.path.join(pwd, args.poscar))
prefix = formulas(structure)[0]
analyzer = SpacegroupAnalyzer(structure, symprec=args.tol)
print(f'Found {prefix}, which has {analyzer.get_space_group_symbol()} '
      f'({str(analyzer.get_space_group_number())}) '
      f'symmetry (with {args.tol} tolerance).'.format())
if args.primitive:
    analyzer = SpacegroupAnalyzer(structure, symprec=0.2)
    primitive = analyzer.find_primitive()
    if len(primitive.sites) != len(structure.sites):
        print('Found primitive structure:')
        print(primitive)
        print('which is now will be used.')
        structure = primitive

if len(args.press) == 1:
    print('Only one pressure is given, output files will be created right here.')
    pressure = arr_str(args.press)
    create_input_opt(args.tol, pwd, structure, args.note, args.o, pressure, args.kppa, args.dyn)
elif len(args.press) > 1:
    print('Multiple pressures are given, creating folder for each pressure...')
    for p in args.press:
        tmp_path = os.path.join(pwd, p)
        if not os.path.isdir(tmp_path):
            os.mkdir(tmp_path)
        note = '_' + p + args.note
        pressure = arr_str(p)
        create_input_opt(args.tol, tmp_path, structure, note, args.o, pressure, args.kppa, args.dyn)

