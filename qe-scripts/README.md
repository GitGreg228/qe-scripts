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