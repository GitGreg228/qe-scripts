from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.cif import CifWriter
from pymatgen.io.vasp.inputs import Kpoints
import numpy as np
import shutil
import json
import os
import re


def parse_mesh(mesh):
    mesh_lst = list()
    res = re.split('(\d+)', mesh)
    for each in res:
        if each.isnumeric():
            _q = int(each)
            mesh_lst.append(_q)
    assert len(mesh_lst) == 3
    return mesh_lst


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


def get_qe_struc(structure, tol, kppa, primitive=True):
    qe_struc = dict()
    qe_struc['prefix'] = formulas(structure)[0]
    analyzer = SpacegroupAnalyzer(structure, symprec=tol)
    sg = analyzer.get_space_group_number()
    qe_struc['sg'] = sg
    qe_struc['sg_str'] = str(sg)
    qe_struc['sym'] = f'{analyzer.get_space_group_symbol()} ({str(sg)})'.format()

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

    if primitive:
        ref_analyzer = SpacegroupAnalyzer(refined, symprec=tol)
        symm_struc = ref_analyzer.get_symmetrized_structure()
        qe_struc['nat'] = str(len(symm_struc.equivalent_sites))
        sites = symm_struc.equivalent_sites
    else:
        ref_analyzer = analyzer
        symm_struc = structure
        qe_struc['nat'] = str(len(symm_struc.sites))
        sites = symm_struc.sites

    qe_struc['ntyp'] = str(len(set(structure.species)))

    species = list()

    for specie in list(set(structure.species)):
        symbol = str(specie.symbol)
        mass = str(specie.atomic_mass).replace(' amu', '')
        upf = symbol + '.UPF'
        species.append(symbol + ' ' + mass + ' ' + upf)
    qe_struc['species'] = "\n".join(species)

    coords = list()

    for site in sites:
        if primitive:
            specie = str(site[0].species)
        else:
            specie = str(site.species)
        atom = ''.join([i for i in specie if not i.isdigit()])
        if primitive:
            x, y, z = site[0].frac_coords
        else:
            x, y, z = site.frac_coords
        coords.append(atom + "\t" + "{:.10f}".format(x) + "\t" + "{:.10f}".format(y) + "\t" + "{:.10f}".format(z))
    qe_struc['positions'] = "\n".join(coords)

    _kpoints = Kpoints.automatic_density_by_vol(refined, kppvol=kppa)

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


def get_contcar(path, o):
    idx = dict()
    idx['s'] = int()
    with open(path, 'r') as f:
        lines = f.readlines()
        f.close()
    for i, line in enumerate(lines):
        if 'Begin final coordinates' in line:
            idx['s'] = i
    for i in range(idx['s'], len(lines)):
        if 'CELL_PARAMETERS' in lines[i]:
            idx['c'] = i
        if 'new unit-cell volume' in lines[i]:
            idx['v'] = i
        if 'g/cm^3' in lines[i] and 'density' in lines[i]:
            idx['d'] = i
        if 'ATOMIC_POSITIONS' in lines[i]:
            idx['a'] = i
        if 'End final coordinates' in lines[i]:
            idx['f'] = i
        if '!    total energy' in lines[i]:
            idx['t'] = i
    if idx['s'] == 0:
        message = f'Warning! No final coordinates found in {path}! Creating CONTCAR from last relaxation step.'.format()
        print(message)
    struc_summ = dict()
    k_Ry_eV = 13.605698066  # 1 Ry = 13.605698066 eV
    nat = idx['f'] - idx['a'] - 1
    vol_A = round(float(lines[idx['v']].split()[-3]), 3)
    E_Ry = round(float(lines[idx['t']].split()[-2]), 3)
    struc_summ['volume'] = [f'{str(vol_A)} Ang^3'.format(),
                            f'{str(round(vol_A / nat, 3))} Ang^3 per atom'.format()]
    struc_summ['density'] = str(round(float(lines[idx['d']].split()[-2]), 3)) + ' ' + lines[idx['d']].split()[-1]
    struc_summ['energy'] = [f'{str(E_Ry)} Ry',
                            f'{str(round(E_Ry / nat, 3))} Ry per atom',
                            f'{str(round(E_Ry * k_Ry_eV, 3))} eV',
                            f'{str(round(E_Ry * k_Ry_eV / nat, 3))} eV per atom']
    # Getting lattice constant
    if 'alat' in lines[idx['c']]:
        alat = lines[idx['c']].split()[-1].replace(')', '')
    else:
        alat = '1'
    # Getting lattice vectors
    lattice = lines[idx['c']+1:idx['a']-1]
    lattice_str = list()
    for vec in lattice:
        _vec = vec.split()
        lattice_str.append('\t'.join(['', *_vec]))
    lattice_str = '\n'.join(lattice_str)
    # Getting atoms positions
    atom_pos = lines[idx['a'] + 1:idx['f']]
    atoms_str = list()
    species = list()
    for a in atom_pos:
        _a = a.split()
        species.append(_a[0])
        atoms_str.append('\t'.join(['', *_a[1:]]))
    atoms_str = '\n'.join(atoms_str)
    # Getting species string
    species_names = list()
    species_numbers = list()
    sp_dict = {i: species.count(i) for i in species}
    prefix = str()
    for sp in sp_dict:
        prefix = prefix + sp
        species_names.append(sp)
        species_numbers.append(str(sp_dict[sp]))
        if sp_dict[sp] != 1:
            prefix = prefix + str(sp_dict[sp])
    volume = ' '.join(struc_summ['volume'])
    density = struc_summ['density']
    energy = ' '.join(struc_summ['energy'])
    prefix_str = f'{prefix} relaxed by QE, V = {volume}, rho = {density}, E = {energy}'

    contcar = [prefix_str, alat, lattice_str,
               '\t'.join(['', *species_names]), '\t'.join(['', *species_numbers]),
               'Direct', atoms_str]
    content = '\n'.join(contcar)

    overwrite(os.path.dirname(path), 'CONTCAR', content, o)

    return struc_summ


def write_json(path, dictionary, name):
    with open(os.path.join(path, name), 'w', encoding='utf-8') as f:
        json.dump(dictionary, f, ensure_ascii=False, sort_keys=True, indent=4)


def print_output(summary):
    message = list()
    for key in summary['symmetry']:
        sym = summary['symmetry'][key]
        message.append(f'{sym} for {key} tolerance')
    message = ', '.join(message)
    print(f'Relaxed structure has the following symmetry: {message}.'.format())
    e = summary['energy'][-3]
    v = summary['volume'][-1]
    d = summary['density']
    print(f'E = {e}, V = {v}, rho = {d}\n'.format())


def reverse_summary(summary):
    reversed_summary = dict()
    for key in summary:
        for _key in summary[key]:
            reversed_summary[_key] = dict()
    for key in summary:
        for _key in summary[key]:
            reversed_summary[_key][key] = summary[key][_key]
    return reversed_summary

