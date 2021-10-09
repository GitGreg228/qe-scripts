from pymatgen.core.structure import IStructure
import argparse
import json
import os

from inputs import *


parser = argparse.ArgumentParser()
parser.add_argument('--path', type=str, default='.', help='Path to a folder with POSCAR')
parser.add_argument('--tol', type=float, default='0.2', help='tolerance')
parser.add_argument('--note', type=str, default='', help='note to add into script')
parser.add_argument('--o', type=boolean_string, default=True, help='Overwrite all files')
parser.add_argument('--tol_max', type=float, default=0.5, help='Maximum allowed tolerance during analysis (NOT CREATING FILES)')
parser.add_argument('--tol_step', type=float, default=0.01, help='Tolerance step during analysis (NOT CREATING FILES)')
parser.add_argument('--kppa', type=int, default=50000, help='K-points per unit volume')
parser.add_argument('--q', nargs='+', default=[], help='Desired q-points meshes for phonon calculations')
parser.add_argument('--save_cif', type=boolean_string, default=True, help='create CIF file')
args = parser.parse_args()

if args.path == '.':
    pwd = str(os.getcwd())
else:
    pwd = os.path.abspath(args.path)

multiple = bool()
tols = dict()

for fname in os.listdir(pwd):
    if 'output.opt' in fname:
        print(f'Found {fname} in {os.path.basename(pwd)}, creating CONTCAR from here')
        get_contcar(os.path.join(pwd, fname), args.o)
        structure = IStructure.from_file(os.path.join(pwd, 'CONTCAR'))
        tols = analyze_symmetry(structure, args.tol_max, args.tol_step, args.save_cif, pwd)
        if len(args.q) > 0:
            create_meshes(args.q, args.tol, pwd, structure,
                          args.note, args.o, args.kppa)
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
                    get_contcar(os.path.join(tmp_path, _fname), args.o)
                    structure = IStructure.from_file(os.path.join(tmp_path, 'CONTCAR'))
                    tols[fname] = analyze_symmetry(structure, args.tol_max, args.tol_step, args.save_cif, tmp_path)
                    if len(args.q) > 0:
                        create_meshes(args.q, args.tol, tmp_path, structure,
                                      args.note, args.o, args.kppa)
    with open(os.path.join(pwd, 'symm.json'), 'w', encoding='utf-8') as f:
        json.dump(tols, f, ensure_ascii=False, indent=4)
        f.close()
