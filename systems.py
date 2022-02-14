import json
import os


def parse_system(path):
    default_system_1 = {
        "partition": "lenovo",
        "modules": "module load intel/mkl-11.2.3 mpi/impi-2018.2.199\nexport PATH=~/qe-6.6/bin/:$PATH",
        "N_pw": 1,
        "n_pw": 8,
        "N_ph": 1,
        "n_ph": 40,
        "pp_path": "../PP",
        "mpirun": "$(which mpirun)"
    }

    default_system_2 = {
        "partition": "all",
        "modules": "module load intel/2020u2\nexport PATH=~/q-e-qe-6.8/bin/:$PATH",
        "N_pw": 1,
        "n_pw": 36,
        "N_ph": 1,
        "n_ph": 36,
        "pp_path": "../PP",
        "mpirun": "srun"
    }

    default_system_3 = {
        "partition": "cpu",
        "modules": "module load compilers/intel-2020\nexport PATH=~/programs/qe-7.0/bin/:$PATH",
        "N_pw": 1,
        "n_pw": 24,
        "N_ph": 2,
        "n_ph": 48,
        "pp_path": "../PP",
        "mpirun": "$(which mpirun)"
    }

    if os.path.isfile(os.path.join(path, 'system.json')):
        with open(os.path.join(path, 'system.json'), 'r') as json_file:
            system = json.load(json_file)
            json_file.close()
    else:
        system = default_system_2
    return system
