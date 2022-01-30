import os


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
    if qe_struc['sg_str'] == '1':
        input_scf = f"""&control
    calculation='scf'
    prefix='{qe_struc['prefix']}',
    pseudo_dir = '.', 
    outdir='.',
 /
&system
    nosym=.TRUE.,
    ibrav= 0, 
    celldm(1) = 0, 
    nat={qe_struc['nat']},
    ntyp={qe_struc['ntyp']},
    ecutwfc=80.0,
    occupations = 'tetrahedra_opt',
 /
 &electrons
 /
ATOMIC_SPECIES
{qe_struc['species']}
ATOMIC_POSITIONS (crystal)
{qe_struc['positions']}
K_POINTS {automatic}
 {kpoints}  {shift}
CELL_PARAMETERS (angstrom)
{qe_struc['lattice']}""".format(automatic=automatic)
    else:
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


def write_ph_in(formula, masses, mesh_lst, q_s, q_f, path, o):
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
  start_q={q_s + 1},
  last_q={q_f + 1},
  nq1 = {mesh_lst[0]}, 
  nq2 = {mesh_lst[1]},
  nq3 = {mesh_lst[2]},
/
""".format()
    if q_s == q_f:
        name = f'ph{q_s + 1}.in'.format()
    else:
        name = f'ph{q_s + 1}to{q_f + 1}.in'.format()
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


def make_3(system, prefix, short, q_s, q_f, len_qpoints, path, o):
    if q_f < len_qpoints - 1:
        _next = f'3{str(q_f + 2)}'.format()
    else:
        _next = '4'
    if q_s == q_f:
        script3 = f"""#!/bin/sh
#SBATCH -o qe.out3 -e qe.err3
#SBATCH -p {system['partition']}
#SBATCH -J {q_s + 1}p{short}
#SBATCH -N {system['N_ph']}
#SBATCH -n {system['n_ph']}

{system['modules']}

echo "PH{q_s + 1} of {prefix} LAUNCHED at" $(date) | tee -a log.{prefix}

{system['mpirun']} $(which ph.x) -in $PWD/ph{q_s + 1}.in &> $PWD/output.ph{q_s + 1}.{prefix} 

echo "PH{q_s + 1} of {prefix} STOPPED at" $(date) | tee -a log.{prefix} 

sbatch script{_next}.sh
""".format()
    else:
        script3 = f"""#!/bin/sh
#SBATCH -o qe.out3 -e qe.err3
#SBATCH -p {system['partition']}
#SBATCH -J {q_s + 1}_{q_f + 1}p{short}
#SBATCH -N {system['N_ph']}
#SBATCH -n {system['n_ph']}

{system['modules']}

echo "PH{q_s + 1} to PH{q_f + 1} of {prefix} LAUNCHED at" $(date) | tee -a log.{prefix}

{system['mpirun']} $(which ph.x) -in $PWD/ph{q_s + 1}to{q_f + 1}.in &> $PWD/output.ph{q_s + 1}to{q_f + 1}.{prefix} 

echo "PH{q_s + 1} to PH{q_f + 1} of {prefix} STOPPED at" $(date) | tee -a log.{prefix} 

sbatch script{_next}.sh
""".format()

    name = f'script3{q_s + 1}.sh'.format()
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

