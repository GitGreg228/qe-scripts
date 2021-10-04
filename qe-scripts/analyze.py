from pymatgen.core.structure import IStructure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.pwscf import PWInput
import numpy as np
import argparse
import shutil
import json
import os
import re

parser = argparse.ArgumentParser()
parser.add_argument('--path', type=str, default='.', help='Path to a folder with POSCAR')
parser.add_argument('--tol_max', type=float, default=0.5, help='Maximum allowed tolerance during analysis (NOT CREATING FILES)')
parser.add_argument('--tol_step', type=float, default=0.01, help='Tolerance step during analysis (NOT CREATING FILES)')
args = parser.parse_args()

if args.path == '.':
    pwd = str(os.getcwd())
else:
    pwd = os.path.abspath(args.path)


def get_contcar(path):
    idx_s = int()
    idx_v = int()
    idx_d = int()
    idx_c = int()
    idx_a = int()
    idx_f = int()
    with open(path, 'r') as f:
        lines = f.readlines()
        f.close()
    for i, line in enumerate(lines):
        if 'Begin final coordinates' in line:
            idx_s = i
        if 0 < idx_s < i and 'CELL_PARAMETERS' in line:
            idx_c = i
        if 0 < idx_s < i and 'new unit-cell volume' in line:
            idx_v = i
        if 0 < idx_s < i and 'g/cm^3' in line and 'density' in line:
            idx_d = i
        if 0 < idx_s < i and 'ATOMIC_POSITIONS' in line:
            idx_a = i
        if 0 < idx_s < i and 'End final coordinates' in line:
            idx_f = i
    if 'alat' in lines[idx_c]:
        alat = lines[idx_c].split()[-1].replace(')', '')
    else:
        alat = '1'
    volume = lines[idx_v].split()[-3] + ' ' + lines[idx_v].split()[-2]
    density = lines[idx_d].split()[-2] + ' ' + lines[idx_d].split()[-1]
    lattice = lines[idx_c+1:idx_a-1]
    lattice_str = list()
    for vec in lattice:
        _vec = vec.split()
        lattice_str.append('\t' + _vec[0] + '\t' + _vec[1] + '\t' + _vec[2] + '\n')
    atom_pos = lines[idx_a+1:idx_f]
    atoms_str = list()
    species = list()
    for a in atom_pos:
        _a = a.split()
        species.append(_a[0])
        atoms_str.append('\t' + _a[1] + '\t' + _a[2] + '\t' + _a[3]+'\n')
    sp_dict = {i: species.count(i) for i in species}
    prefix = str()
    species_str = ['', '']
    for sp in sp_dict:
        prefix = prefix + sp
        species_str[0] = species_str[0] + '\t' + sp
        species_str[1] = species_str[1] + '\t' + str(sp_dict[sp])
        if sp_dict[sp] != 1:
            prefix = prefix + str(sp_dict[sp])
    prefix_str = f'{prefix} relaxed by QE, V = {volume}, rho = {density}'
    with open(os.path.join(path, '../CONTCAR'), 'w') as f:
        f.write(prefix_str + '\n' + alat + '\n')
        for vec in lattice_str:
            f.write(vec)
        f.write(species_str[0] + '\n' + species_str[1] + '\nDirect\n')
        for a in atoms_str:
            f.write(a)
        f.close()


def analyze_symmetry(structure):
    prev_number = 0
    tols = dict()
    for tol in np.arange(0.01, args.tol_max, args.tol_step):
        analyzer = SpacegroupAnalyzer(structure, symprec=tol)
        number = analyzer.get_space_group_number()
        if number > prev_number:
            symbol = analyzer.get_space_group_symbol()
            tols["{:0.2f}".format(tol)] = str(number) + '(' + symbol + ')'
        prev_number = number
    return tols


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


def create_input_scf(tol, path, structure, note, o, kppa):
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

    input_scf = f"""&control
    calculation='scf'
    prefix='{prefix}',
    pseudo_dir = '.', 
    outdir='.',
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
        occupations = 'tetrahedra_opt',
 /
 &electrons
 /
ATOMIC_SPECIES
{species}
ATOMIC_POSITIONS (crystal_sg)
{positions}
K_POINTS {automatic}
 {kpoints}  {shift}""".format(automatic=automatic)

    overwrite(path, 'input.scf', input_scf, o)

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

    script2 = f"""#!/bin/sh
#SBATCH -o qe.out1 -e qe.err1
#SBATCH -p {system['partition']}
#SBATCH -J s{short}
#SBATCH -N {system['N_pw']}
#SBATCH -n {system['n_pw']}

{system['modules']}

echo "SCF of {prefix} LAUNCHED at" $(date) | tee -a log.{prefix} 

$(which mpirun) $(which pw.x) -in $PWD/input.scf &> $PWD/output.scf.{prefix} 

echo "SCF of {prefix} STOPPED at" $(date) | tee -a log.{prefix} 

""".format()

    overwrite(path, 'script2.sh', script2, o)


multiple = bool()
tols = dict()

for fname in os.listdir(pwd):
    if 'output.opt' in fname:
        print(f'Found {fname} in {os.path.basename(pwd)}, creating CONTCAR from here')
        get_contcar(os.path.join(pwd, fname))
        structure = IStructure.from_file(os.path.join(pwd, 'CONTCAR'))
        tols = analyze_symmetry(structure)
        with open(os.path.join(pwd, 'symm.json'), 'w', encoding='utf-8') as f:
            json.dump(tols, f, ensure_ascii=False, indent=4)
            f.close()
        multiple = False
    else:
        multiple = True

if multiple:
    dir_lst = list()
    for fname in os.listdir(pwd):
        tmp_path = os.path.join(pwd, fname)
        if os.path.isdir(tmp_path):
            for _fname in os.listdir(tmp_path):
                if 'output.opt' in _fname:
                    print(f'Found {_fname} in {os.path.basename(tmp_path)}, creating CONTCAR from here')
                    get_contcar(os.path.join(tmp_path, _fname))
                    structure = IStructure.from_file(os.path.join(tmp_path, 'CONTCAR'))
                    tols[fname] = analyze_symmetry(structure)
    with open(os.path.join(pwd, 'symm.json'), 'w', encoding='utf-8') as f:
        json.dump(tols, f, ensure_ascii=False, indent=4)
        f.close()