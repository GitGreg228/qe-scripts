from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.cif import CifWriter
from pymatgen.io.vasp.inputs import Kpoints
import numpy as np
import shutil
import os
import re


def boolean_string(s):
    if s not in {'False', 'True'}:
        raise ValueError('Not a valid boolean string')
    return s == 'True'


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


def arr_str(s):
    return str(s).replace('[', '').replace(']', '').replace(',', '').replace('\'', '')


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


def formulas(structure):
    formula = structure.formula.replace(' ', '')
    res = re.split('(\d+)', formula)
    prefix = str()
    short = str()
    for c in res:
        if c.isalpha():
            prefix = prefix + c
        elif c.isnumeric():
            short = short + c
            if not c == '1':
                prefix = prefix + c
    return prefix, short


def get_masses(structure):
    masses = str()
    species = list(set(structure.species))
    for i in range(len(species)):
        mass = str(species[i].atomic_mass).replace(' amu', '')
        masses = masses + f'amass({i+1})={mass},\n  '.format()
    return masses


def analyze_symmetry(structure, tol_max, tol_step, save_cif, path):
    prev_number = 0
    tols = dict()
    for tol in np.arange(0.01, tol_max, tol_step):
        analyzer = SpacegroupAnalyzer(structure, symprec=tol)
        number = analyzer.get_space_group_number()
        if number > prev_number:
            symbol = analyzer.get_space_group_symbol()
            tols["{:0.2f}".format(tol)] = str(number) + '(' + symbol + ')'
            if save_cif:
                name = f'{formulas(structure)[0]}_{os.path.basename(path)}' \
                       f'_{analyzer.get_space_group_number()}.cif'.format()
                CifWriter(analyzer.get_symmetrized_structure(), symprec=tol).write_file(os.path.join(path, name))
        prev_number = number
    return tols


def get_qe_struc(structure, tol, kppa):
    qe_struc = dict()
    qe_struc['prefix'] = formulas(structure)[0]
    analyzer = SpacegroupAnalyzer(structure, symprec=tol)
    sg = analyzer.get_space_group_number()
    qe_struc['sg_str'] = str(sg)

    if sg < 142:
        qe_struc['uniqueb'] = '\nuniqueb = .TRUE.,'
        refined = analyzer.get_refined_structure()
    else:
        qe_struc['uniqueb'] = ''
        refined = analyzer.get_primitive_standard_structure()
    qe_struc['refined'] = refined
    
    qe_struc['a'] = round(refined.lattice.a, 10)
    qe_struc['b'] = round(refined.lattice.b, 10)
    qe_struc['c'] = round(refined.lattice.c, 10)
    qe_struc['cosBC'] = round(np.cos(refined.lattice.alpha * np.pi / 180), 10)
    qe_struc['cosAC'] = round(np.cos(refined.lattice.beta * np.pi / 180), 10)
    qe_struc['cosAB'] = round(np.cos(refined.lattice.gamma * np.pi / 180), 10)

    ref_analyzer = SpacegroupAnalyzer(refined, symprec=tol)
    symm_struc = ref_analyzer.get_symmetrized_structure()

    qe_struc['nat'] = str(len(symm_struc.equivalent_sites))
    qe_struc['ntyp'] = str(len(set(structure.species)))

    species = list()

    for specie in list(set(structure.species)):
        symbol = str(specie.symbol)
        mass = str(specie.atomic_mass).replace(' amu', '')
        upf = symbol + '.UPF'
        species.append(symbol + ' ' + mass + ' ' + upf)
    qe_struc['species'] = "\n".join(species)

    coords = list()

    eq_sites = symm_struc.equivalent_sites
    for site in eq_sites:
        specie = str(site[0].species)
        atom = ''.join([i for i in specie if not i.isdigit()])
        x = site[0].frac_coords[0]
        y = site[0].frac_coords[1]
        z = site[0].frac_coords[2]
        coords.append(atom + "\t" + "{:.10f}".format(x) + "\t" + "{:.10f}".format(y) + "\t" + "{:.10f}".format(z))
    qe_struc['positions'] = "\n".join(coords)

    _kpoints = Kpoints.automatic_density(refined, kppa=kppa)

    qe_struc['kpoints'] = arr_str(_kpoints.kpts[0])
    qe_struc['shift'] = arr_str(_kpoints.kpts_shift)
    return qe_struc


def copy_pp(system, path):
    if os.path.isdir(system['pp_path']):
        for fname in os.listdir(system['pp_path']):
            if '.UPF' in fname:
                tmp_src = os.path.join(system['pp_path'], fname)
                if os.path.isfile(tmp_src):
                    shutil.copy2(tmp_src, path)


def parse_mesh(mesh):
    mesh_lst = list()
    res = re.split('(\d+)', mesh)
    for each in res:
        if each.isnumeric():
            _q = int(each)
            mesh_lst.append(_q)
    assert len(mesh_lst) == 3
    return mesh_lst
