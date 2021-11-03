from pymatgen.core.structure import IStructure
import argparse

from inputs import *

"""
Parsing arguments
"""

parser = argparse.ArgumentParser()
parser.add_argument('--tol', type=float, default='0.2', help='tolerance')
parser.add_argument('--path', type=str, default='.', help='Path to a folder with POSCAR')
parser.add_argument('--poscar', type=str, default='POSCAR', help='POSCAR file name')
parser.add_argument('--note', type=str, default='', help='note to add into script')
parser.add_argument('--o', type=boolean_string, default=True, help='Overwrite all files')
parser.add_argument('--press', nargs='+', default=[2000], help='Pressure(s) in kBar')
parser.add_argument('--kppa', type=int, default=300, help='K-points per unit volume')
parser.add_argument('--primitive', default=True, type=boolean_string, help='if print primitive structure')
parser.add_argument('--dyn', type=str, default='all', help='cell dynamics')
args = parser.parse_args()

if args.path == '.':
    pwd = str(os.getcwd())
else:
    pwd = os.path.abspath(args.path)


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
        note = '_' + p + '_' + args.note
        pressure = arr_str(p)
        create_input_opt(args.tol, tmp_path, structure, note, args.o, pressure, args.kppa, args.dyn, args.primitive)

