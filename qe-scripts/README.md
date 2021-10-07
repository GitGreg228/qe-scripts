# Quantum Espresso scripts

## `input.opt` generator

Currently, several QE input generators are available on the Internet (for example, https://www.materialscloud.org/work/tools/qeinputgenerator). However, for some reason, they are not capable in detecting crystal symmetry, which can crucially affect the calculation speed.
`qe_sc.py` is used to set up the `vc-relax` calculation for given structure. It generates `input.opt` and `script1.sh` files basing on given POSCAR and (optionally) system.json file. 
* `input.opt` contains parameters of calculation and crystal structure geometry. Only inequvalently placed atoms are written in this file. The `space_group` parameter defines the symmetry, and the `A`, `B`, `C`, `cosAB`, `cosAC`, `cosBC` define the unit cell shape and volume. 
* `script1.sh` is a script for Slurm Workload Manager

### Usage

For `input.opt` generation, run the 
```bash
python qe_sc.py
```
in a folder containing POSCAR

User can choose the tolerance running 
```bash
python qe_sc.py --tol=%desired_tolerance%
```
Default tolerance is 0.2.

The pressure can also be set up by
```bash
python qe_sc.py --press=%desired_pressure_in_kbar%
```
Default pressure is 2000.

If the POSCAR placed in another path (for example, in `~/QE/H3S`) and named differently (for example, `CONTCAR`), user can run
```bash
python qe_sc.py --path=~/QE/H3S --poscar=CONTCAR
```

### Multiple pressures

The program `qe_sc.py` allows to create multiple folders for each given pressure. To set up multiple calculations, run
```bash
python qe_sc.py --press %P1% %P2% %P3% ... %Pn% 
```
As result, corresponding folders will appear in the folder next to `POSCAR` file.

## Symmetry analyzer

Symmetry analyzer `analyze.py` code allows to conveniently analyze symmetry of relaxed structures. It looks for `output.opt.` files next to `POSCAR` (or in subdirectories next to `POSCAR`, which could be created by `qe_sc.py`). When `output.opt.` is found, final lattice and atomic coordinates are parsed. Then the `CONTCAR` and the symmetrized `.cif` file are generated in the folder. The summary of all found symmetries are placed in `symm.json` file in the directory in which `analyze.py` was launched.
The algorithm of symmetry analysis is the same as in [`symmetrize.py`](https://github.com/GitGreg228/cms-scripts/tree/main/vasp-scripts)

### Usage

In a directory containing `output.opt.` or containing subdirectories with `output.opt.`, run
```bash
python analyze.py
```
You can choose range of tolerances in which the structure will be analyzed.
```bash
python analyze.py --tol_step=%desired_tolerance_step% --tol_max=%desired_max_tolerance%
```
Default `--tol_step` and `--tol_max` are 0.01 and 0.5 accordingly.
You can disable creation of `.cif` files using
