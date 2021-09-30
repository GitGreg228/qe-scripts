# VASP scripts

## Crystal structure symmetrizer

Common crystal structure symmetrizers (for example, https://uspex-team.org/online_utilities/poscar2cif/) allow to symmetrize a crystal with specified tolerances, however, they are inconvenient if the task is to see how the symmetry changes with tolerance parameter.

The script `symmetrize.py` analyzes the symmetry on given range of tolerances with a given step.

### How to use

Run
```bash
python symmetrize.py
```
in a folder containing POSCAR file. As result, the `symm.json` file will be created, containing a dictionary with tolerances and space groups.

### Options

If your POSCAR file is named differently (for example, CONTCAR), use:
```bash
python symmetrize.py --poscar=CONTCAR
```
To choose another folder, run
```bash
python symmetrize.py --path=your/desired/path
```
to change the range and step of tolerances, use `--tol_step` and `--tol_max` parameters:
```bash
python symmetrize.py --tol_step=0.05 --tol_max=0.2
```
Default `--tol_step` and `--tol_max` are 0.01 and 0.5 accordingly.
You can also create folders with symmetrized POSCARS and cifs using
```bash
python symmetrize.py --q=False
```