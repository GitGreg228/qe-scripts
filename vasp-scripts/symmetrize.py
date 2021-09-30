from pymatgen.core.structure import IStructure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.cif import CifWriter, CifParser
import numpy as np
import argparse
import json
import os

parser = argparse.ArgumentParser()
parser.add_argument('--tol_max', type=float, default=0.5, help='Maximum allowed tolerance')
parser.add_argument('--tol_step', type=float, default=0.01, help='Tolerance step')
parser.add_argument('--path', type=str, default='.', help='Path to a folder with POSCAR')
parser.add_argument('--poscar', type=str, default='POSCAR', help='POSCAR file name')
parser.add_argument('--q', type=bool, default=False, help='Quiet mode, create only symm.json')

args = parser.parse_args()

if args.path == '.':
    pwd = str(os.getcwd())
else:
    pwd = os.path.abspath(args.path)

structure = IStructure.from_file(os.path.join(pwd, args.poscar))
prev_number = 0
tols = dict()

if args.q:
    dir_name = os.path.basename(pwd) + '_0'
    folder = os.path.join(pwd, dir_name)
    if not os.path.exists(folder):
        os.mkdir(folder)
    Poscar(structure).write_file(os.path.join(folder, args.poscar))

for tol in np.arange(0.01, args.tol_max, args.tol_step):
    analyzer = SpacegroupAnalyzer(structure, symprec=tol)
    number = analyzer.get_space_group_number()
    if number > prev_number:
        symbol = analyzer.get_space_group_symbol()
        lat = analyzer.get_lattice_type()
        tols["{:0.2f}".format(tol)] = str(number) + '(' + symbol + ')'
        if args.q:
            dir_name = os.path.basename(pwd) + '_' + str(number) + '_' + lat
            folder = os.path.join(pwd, dir_name)
            if not os.path.exists(folder):
              os.mkdir(folder)
            symm_struc = analyzer.get_symmetrized_structure()
            ref_struc = analyzer.get_refined_structure()
            Poscar(ref_struc).write_file(os.path.join(folder, args.poscar))
            name = symm_struc.formula.replace(' ', '') + '_' + str(number) + '.cif'
            path = os.path.join(folder, name)
            CifWriter(symm_struc, symprec=tol).write_file(path)
            cif_struc = CifParser(path).get_structures()[0]
    prev_number = number

with open(os.path.join(pwd, 'symm.json'), 'w', encoding='utf-8') as f:
    json.dump(tols, f, ensure_ascii=False, indent=4)
