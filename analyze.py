from pymatgen.core.structure import IStructure
import argparse

from utils import *


parser = argparse.ArgumentParser()
parser.add_argument('--path', type=str, default='.', help='Path to a folder with POSCAR')
parser.add_argument('--tol', type=float, default='0.2', help='tolerance')
parser.add_argument('--note', type=str, default='', help='note to add into script')
parser.add_argument('--o', type=boolean_string, default=True, help='Overwrite all files')
parser.add_argument('--tol_max', type=float, default=0.5, help='Maximum allowed tolerance during analysis (NOT CREATING FILES)')
parser.add_argument('--tol_step', type=float, default=0.01, help='Tolerance step during analysis (NOT CREATING FILES)')
parser.add_argument('--kppa', type=int, default=2000, help='Grid density per Angstrom^(-3)')
parser.add_argument('--q', nargs='+', default=[], help='Desired q-points meshes for phonon calculations')
parser.add_argument('--q_num', type=int, default=False, help='Desired q-points meshes for phonon calculations')
parser.add_argument('--save_cif', type=boolean_string, default=True, help='create CIF file')
parser.add_argument('--mul', type=int, default=4, help='k = q * multiplier')
parser.add_argument('--subs', type=boolean_string, default=False, help='Launch multiple subsequent scripts in ph calculation')
parser.add_argument('--nosym', type=boolean_string, default=True, help='Use no symmetry after relaxation')
args = parser.parse_args()

if args.path == '.':
    pwd = str(os.getcwd())
else:
    pwd = os.path.abspath(args.path)

multiple = bool()
summary = dict()
qpoints = dict()

for fname in os.listdir(pwd):
    if 'output.opt' in fname:
        print(f'Found {fname} in {os.path.basename(pwd)}, creating CONTCAR from here')
        summary = get_contcar(os.path.join(pwd, fname), args.o)
        structure = IStructure.from_file(os.path.join(pwd, 'CONTCAR'))
        summary['symmetry'] = analyze_symmetry(structure, args.tol_max, args.tol_step, args.save_cif, pwd)
        print_output(summary)
        if len(args.q) > 0:
            qpoints = create_meshes(args.q, args.tol, pwd, structure,
                                    args.note, args.o, args.kppa, args.mul, args.subs, args.nosym,  q_num=args.q_num)
            write_json(pwd, qpoints, 'qpoints.json')
        write_json(pwd, summary, 'summary.json')
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
                    print(f'Found {_fname} in {os.path.basename(tmp_path)}, creating CONTCAR from here...')
                    summary[fname] = get_contcar(os.path.join(tmp_path, _fname), args.o)
                    structure = IStructure.from_file(os.path.join(tmp_path, 'CONTCAR'))
                    summary[fname]['symmetry'] = analyze_symmetry(structure, args.tol_max,
                                                                  args.tol_step, args.save_cif, tmp_path)
                    print_output(summary[fname])
                    if len(args.q) > 0:
                        qpoints[fname] = create_meshes(args.q, args.tol, tmp_path, structure,
                                                       args.note, args.o, args.kppa, args.mul,
                                                       args.subs, args.nosym, q_num=args.q_num)
    if len(args.q) > 0:
        write_json(pwd, qpoints, 'qpoints.json')
    write_json(pwd, reverse_summary(summary), 'summary.json')

print('Parameters of relaxed structure(s) are saved in summary.json')
if len(args.q) > 0:
    print('Parameters of generated qpoints for phonon calculations are saved in qpoints.json')
