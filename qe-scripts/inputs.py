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


def write_ph_in(formula, masses, mesh_lst, q, path, o):
    ph_in = f"""Electron-phonon coefficients for {formula}
 &inputph
  tr2_ph=1.0d-8,
  prefix='{formula}',
  outdir = '.',
  fildvscf='dv',
  fildyn = '{formula}.dyn',
  {masses}fildrho = 'drho',
  ldisp = .true.,
  lshift_q = .true.,
  start_q={q + 1},
  last_q={q + 1},
  nq1 = {mesh_lst[0]}, 
  nq2 = {mesh_lst[1]},
  nq3 = {mesh_lst[2]},
/
""".format()

    name = f'ph{q + 1}.in'.format()
    overwrite(path, name, ph_in, o)


def write_elph_in(masses, mesh_lst, len_qpoints, kpoints, prefix, path, o):
    elph_in = f""" Electron-phonon for {prefix}
&INPUTPH
tr2_ph=1.0d-8,
prefix = '{prefix}',
outdir = '.',
fildvscf = 'dv',
fildyn = '{prefix}.dyn'
{masses}fildrho = 'drho',
ldisp = .true.,
lshift_q = .true.,
start_q=1,
last_q={len_qpoints},
nq1 = {mesh_lst[0]}, 
nq2 = {mesh_lst[1]},
nq3 = {mesh_lst[2]},
electron_phonon = "lambda_tetra"
nk1 = {kpoints.split()[0]},
nk2 = {kpoints.split()[1]},
nk3 = {kpoints.split()[2]},
/
&INPUTa2F
nfreq = 10000
/""".format()
    
    overwrite(path, 'elph.in', elph_in, o)


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

sbatch script31.sh
""".format()

    overwrite(path, 'script2.sh', script2, o)


def make_3(system, prefix, short, q, len_qpoints, path, o):
    if q < len_qpoints - 1:
        _next = f'3{str(q + 2)}'.format()
    else:
        _next = '4'
    script3 = f"""#!/bin/sh
#SBATCH -o qe.out3 -e qe.err3
#SBATCH -p {system['partition']}
#SBATCH -J {q + 1}p{short}
#SBATCH -N {system['N_ph']}
#SBATCH -n {system['n_ph']}

{system['modules']}

echo "PH{q + 1} of {prefix} LAUNCHED at" $(date) | tee -a log.{prefix}

{system['mpirun']} $(which ph.x) -in $PWD/ph{q + 1}.in &> $PWD/output.ph{q + 1}.{prefix} 

echo "PH{q + 1} of {prefix} STOPPED at" $(date) | tee -a log.{prefix} 

sbatch script{_next}.sh
""".format()

    name = f'script3{q + 1}.sh'.format()
    overwrite(path, name, script3, o)


def make_4(system, prefix, short, path, o):
    script4 = f"""#!/bin/sh
#SBATCH -o qe.out4 -e qe.err4
#SBATCH -p {system['partition']}
#SBATCH -J e{short}
#SBATCH -N {system['N_ph']}
#SBATCH -n {system['n_ph']}

{system['modules']}

echo "ELPH of {prefix} LAUNCHED at" $(date) | tee -a log.{prefix}

{system['mpirun']} $(which ph.x) -in $PWD/elph.in &> $PWD/output.elph.{prefix} 

echo "ELPH of {prefix} STOPPED at" $(date) | tee -a log.{prefix} 

""".format()

    overwrite(path, 'script4.sh', script4, o)


def create_input_opt(tol, path, structure, note, o, pressure, kppa, dyn):
    qe_struc = get_qe_struc(structure, tol, kppa)
    write_opt(qe_struc, pressure, dyn, path, o)
    system = parse_system(path)
    copy_pp(system, path)
    prefix = formulas(structure)[0]
    short = formulas(structure)[1] + f'{note}'.format()
    make_1(system, prefix, short, path, o)


def create_input_scf(path, qe_struc, o, kpoints, prefix, short, shift='0 0 0'):
    write_scf(qe_struc, kpoints, shift, path, o)
    system = parse_system(path)
    copy_pp(system, path)
    make_2(system, prefix, short, path, o)


def create_ph_ins(path, mesh, mesh_lst, structure, kpoints, tol, prefix, short, o):
    analyzer = SpacegroupAnalyzer(structure, symprec=tol)
    qpoints = analyzer.get_ir_reciprocal_mesh(tuple(mesh_lst), is_shift=(0.5, 0.5, 0.5))
    masses = get_masses(structure)
    len_qpoints = len(qpoints)
    system = parse_system(path)
    for q in range(len_qpoints):
        write_ph_in(prefix, masses, mesh_lst, q, path, o)
        make_3(system, prefix, short, q, len_qpoints, path, o)
    qpoints_dict = dict()
    for qpoint in qpoints:
        xyz = np.round(qpoint[0], 3).tolist()
        qpoints_dict['  '.join(["%.3f" % k for k in xyz])] = qpoint[1].item()
    make_4(system, prefix, short, path, o)
    write_elph_in( masses, mesh_lst, len_qpoints, kpoints, prefix, path, o)
    return qpoints_dict


def create_meshes(q, tol, path, structure, note, o, kppa):
    mesh_dict = dict()
    for mesh in q:
        _tmp_path = os.path.join(path, mesh)
        if not os.path.isdir(_tmp_path):
            os.mkdir(_tmp_path)
        kpoints = list()
        refined = get_qe_struc(structure, tol, kppa)['refined']
        _kpoints = Kpoints.automatic_density_by_vol(refined, kppvol=kppa).kpts[0]
        mesh_lst = parse_mesh(mesh)
        k_total = 1
        q_total = 1
        for i, _q in enumerate(mesh_lst):
            kpoints.append(str(ceil(_kpoints[i]/_q)*_q))
            k_total = k_total * _kpoints[i]
            q_total = q_total * _q
        kpoints = ' '.join(kpoints)
        if note:
            note = f'_{os.path.basename(path)}_{mesh}_{note}'.format()
        else:
            note = f'_{os.path.basename(path)}_{mesh}'.format()
        qe_struc = get_qe_struc(structure, tol, kppa=1)
        prefix = formulas(structure)[0]
        short = formulas(structure)[1] + f'{note}'.format()
        create_input_scf(_tmp_path, qe_struc, o, kpoints, prefix, short)
        mesh_dict[mesh] = dict()
        mesh_dict[mesh]['space_group'] = qe_struc['sym']
        mesh_dict[mesh]['xyz'] = create_ph_ins(_tmp_path, mesh, mesh_lst, structure, kpoints, tol, prefix, short, o)
        mesh_dict[mesh]['kpoints_scf'] = kpoints
        mesh_dict[mesh]['vol_Ang^3'] = round(structure.volume, 3)
        mesh_dict[mesh]['vol_reciprocal_Ang^(-3)'] = k_total / kppa
        mesh_dict[mesh]['kppvol'] = kppa
        mesh_dict[mesh]['qppa'] = round(kppa * q_total / k_total)
        mesh_dict[mesh]['qppu'] = q_total
        mesh_dict[mesh]['qppuirr'] = len(mesh_dict[mesh]['xyz'])

    return mesh_dict
