import json
import os


def parse_system(path):
    default_system_1 = {
        "partition": "lenovo",
        "modules": "module load intel2021/tbb/2021.2.0 intel2021/compiler-rt/2021.2.0 intel2021/mpi/2021.2.0 intel2021/debugger/10.1.1 intel2021/compiler/2021.2.0 intel2021/mkl/2021.2.0\nexport PATH=~/programs/qe-7.0/bin/:$PATH",
        "N_pw": 1,
        "n_pw": 8,
        "N_ph": 1,
        "n_ph": 20,
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
