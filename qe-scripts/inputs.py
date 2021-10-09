from math import ceil
from utils import *

from systems import parse_system


def write_opt(qe_struc, pressure, dyn, path, o):
    """
    Writing input.opt
    """
    automatic = '{automatic}'

    input_opt = f"""&control
    calculation='vc-relax'
    restart_mode='from_scratch',
    prefix={qe_struc['prefix']},
    pseudo_dir = '.',
    outdir='.',
    wf_collect = .true.,
    nstep = 9999,
/
&system
        space_group={qe_struc['sg_str']},
        A={qe_struc['a']},
        B={qe_struc['b']},
        C={qe_struc['c']},
        cosAB={qe_struc['cosAB']},
        cosBC={qe_struc['cosBC']},
        cosAC={qe_struc['cosAC']},{qe_struc['uniqueb']}
        nat={qe_struc['nat']},
        ntyp={qe_struc['ntyp']},
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
{qe_struc['species']}
ATOMIC_POSITIONS (crystal_sg)
{qe_struc['positions']}
K_POINTS {automatic}
 {qe_struc['kpoints']}  {qe_struc['shift']}""".format(automatic=automatic)

    overwrite(path, 'input.opt', input_opt, o)


def write_scf(qe_struc, kpoints, shift, path, o):
    """
    Writing input.scf
    """
    automatic = '{automatic}'

    input_scf = f"""&control
    calculation='scf'
    prefix='{qe_struc['prefix']}',
    pseudo_dir = '.', 
    outdir='.',
 /
&system
        space_group={qe_struc['sg_str']},
        A={qe_struc['a']},
        B={qe_struc['b']},
        C={qe_struc['c']},
        cosAB={qe_struc['cosAB']},
        cosBC={qe_struc['cosBC']},
        cosAC={qe_struc['cosAC']},{qe_struc['uniqueb']}
        nat={qe_struc['nat']},
        ntyp={qe_struc['ntyp']},
        ecutwfc=80.0,
        occupations = 'tetrahedra_opt',
 /
 &electrons
 /
ATOMIC_SPECIES
{qe_struc['species']}
ATOMIC_POSITIONS (crystal_sg)
{qe_struc['positions']}
K_POINTS {automatic}
 {kpoints}  {shift}""".format(automatic=automatic)

    overwrite(path, 'input.scf', input_scf, o)


def make_1(system, prefix, short, path, o):
    script1 = f"""#!/bin/sh
#SBATCH -o qe.out1 -e qe.err1
#SBATCH -p {system['partition']}
#SBATCH -J r{short}
#SBATCH -N {system['N_pw']}
#SBATCH -n {system['n_pw']}

{system['modules']}

echo "OPT of {prefix} LAUNCHED at" $(date) | tee -a log.{prefix} 

{system['mpirun']} $(which pw.x) -in $PWD/input.opt &> $PWD/output.opt.{prefix} 

echo "OPT of {prefix} STOPPED at" $(date) | tee -a log.{prefix} 

""".format()

    overwrite(path, 'script1.sh', script1, o)


def make_2(system, prefix, short, path, o):
    script2 = f"""#!/bin/sh
#SBATCH -o qe.out2 -e qe.err2
#SBATCH -p {system['partition']}
#SBATCH -J s{short}
#SBATCH -N {system['N_pw']}
#SBATCH -n {system['n_pw']}

{system['modules']}

echo "SCF of {prefix} LAUNCHED at" $(date) | tee -a log.{prefix}

{system['mpirun']} $(which pw.x) -in $PWD/input.scf &> $PWD/output.scf.{prefix} 

echo "SCF of {prefix} STOPPED at" $(date) | tee -a log.{prefix} 

""".format()

    overwrite(path, 'script2.sh', script2, o)


def create_input_opt(tol, path, structure, note, o, pressure, kppa, dyn):
    qe_struc = get_qe_struc(structure, tol, kppa)
    write_opt(qe_struc, pressure, dyn, path, o)
    system = parse_system(path)
    copy_pp(system, path)
    prefix = formulas(structure)[0]
    short = formulas(structure)[1] + f'{note}'.format()
    make_1(system, prefix, short, path, o)


def create_input_scf(tol, path, qe_struc, note, o, kpoints, prefix, short, shift='0 0 0'):
    write_scf(qe_struc, kpoints, shift, path, o)
    system = parse_system(path)
    copy_pp(system, path)
    make_2(system, prefix, short, path, o)


def create_meshes(q, tol, path, structure, note, o, multiplier, kppa):
    for mesh in q:
        _tmp_path = os.path.join(path, mesh)
        if not os.path.isdir(_tmp_path):
            os.mkdir(_tmp_path)
        res = re.split('(\d+)', mesh)
        mesh_lst = list()
        for each in res:
            if each.isnumeric():
                mesh_lst.append(each)
        assert len(mesh_lst) == 3
        qpoints = str()
        for _q in mesh_lst:
            qpoints = qpoints + str(multiplier * int(_q)) + ' '
        kpoints = str()
        refined = get_qe_struc(structure, tol, kppa)['refined']
        _kpoints = Kpoints.automatic_density(refined, kppa=kppa).kpts[0]
        for i, _q in enumerate(mesh_lst):
            kpoints = kpoints + str(ceil(_kpoints[i]/int(_q))*int(_q)) + ' '
        note = f'_{os.path.basename(path)}_{mesh}'.format()
        qe_struc = get_qe_struc(structure, tol, kppa=1)
        prefix = formulas(structure)[0]
        short = formulas(structure)[1] + f'{note}'.format()
        create_input_scf(tol, _tmp_path, qe_struc, note, o, kpoints, prefix, short)


def get_contcar(path, o):
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
    for i, line in enumerate(lines, start=idx_s):
        if 'CELL_PARAMETERS' in line:
            idx_c = i
        if 'new unit-cell volume' in line:
            idx_v = i
        if 'g/cm^3' in line and 'density' in line:
            idx_d = i
        if 'ATOMIC_POSITIONS' in line:
            idx_a = i
        if 'End final coordinates' in line:
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
        lattice_str.append('\t' + _vec[0] + '\t' + _vec[1] + '\t' + _vec[2])
    lattice_str = "\n".join(lattice_str)
    atom_pos = lines[idx_a+1:idx_f]
    atoms_str = list()
    species = list()
    for a in atom_pos:
        _a = a.split()
        species.append(_a[0])
        atoms_str.append('\t' + _a[1] + '\t' + _a[2] + '\t' + _a[3])
    atoms_str = "\n".join(atoms_str)
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
    content = prefix_str + '\n' + alat + '\n' + lattice_str + \
              '\n' + species_str[0] + '\n' + species_str[1] + '\nDirect\n' + atoms_str
    overwrite(os.path.dirname(path), 'CONTCAR', content, o)
