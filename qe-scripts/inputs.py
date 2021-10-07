import json
from math import ceil
from utils import *


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


def parse_system(path):
    default_system = {
        "partition": "lenovo",
        "modules": "module load intel/mkl-11.2.3 mpi/impi-2018.2.199\nexport PATH=~/qe-6.6/bin/:s$PATH",
        "N_pw": 1,
        "n_pw": 8,
        "N_ph": 1,
        "n_ph": 40,
        "pp_path": "../../PP"
    }

    if os.path.isfile(os.path.join(path, 'system.json')):
        with open(os.path.join(path, 'system.json'), 'r') as json_file:
            system = json.load(json_file)
            json_file.close()
    else:
        system = default_system
    return system


def make_1(system, prefix, short, path, o):
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


def make_2(system, prefix, short, path, o):
    script2 = f"""#!/bin/sh
#SBATCH -o qe.out2 -e qe.err2
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


def create_input_opt(tol, path, structure, note, o, pressure, kppa, dyn):
    qe_struc = get_qe_struc(structure, tol, kppa)
    write_opt(qe_struc, pressure, dyn, path, o)
    system = parse_system(path)
    copy_pp(system, path)
    prefix = formulas(structure)[0]
    short = formulas(structure)[1] + f'{note}'.format()
    make_1(system, prefix, short, path, o)


def create_input_scf(tol, path, structure, note, o, kpoints, shift='0 0 0'):
    qe_struc = get_qe_struc(structure, tol, kppa=1)
    write_scf(qe_struc, kpoints, shift, path, o)
    system = parse_system(path)
    copy_pp(system, path)
    prefix = formulas(structure)[0]
    short = formulas(structure)[1] + f'{note}'.format()
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
        _kpoints = Kpoints.automatic_density(structure, kppa=kppa).kpts[0]
        for i, _q in enumerate(mesh_lst):
            kpoints = kpoints + str(ceil(_kpoints[i]/int(_q))*int(_q)) + ' '
        create_input_scf(tol, _tmp_path, structure, note, o, kpoints)
