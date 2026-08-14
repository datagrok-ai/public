from rdkit.Chem import Descriptors
from rdkit.Chem import Descriptors3D
from rdkit.Chem import AllChem
import rdkit.Chem as Chem
import numpy as np
import sys

from utils import np_none, set_type

_descriptors_tree = None
_descriptors_list = None
_descriptor_types = {}

_module_names = {'EState.EState_VSA': 'EState VSA', 'EState.EState': 'EState', 'Descriptors3D': 'Descriptors 3D'}


def _is_3d(func):
    return func.__module__ == 'rdkit.Chem.Descriptors3D'


def _first_doc_line(obj):
    for line in (obj.__doc__ or '').splitlines():
        line = line.strip()
        if line and ') -> ' not in line:
            return line
    return ''


def get_3d_descriptors_names():
    """
    Gets names of the descriptors derived from a molecule's 3D structure.

    :return: List of 3D descriptors names.
    """
    return [name for name, func in _get_descriptors_list() if _is_3d(func)]


def get_descriptors_tree():
    """
    Gets all available Molecule descriptors in tree view (module -> descriptor).

    :return: List of available Molecule descriptors.
    """
    global _descriptors_tree
    if _descriptors_tree is None:
        _descriptors_tree = {}
        for name, func in _get_descriptors_list():
            module_name = func.__module__.replace('rdkit.Chem.', '')
            descriptor = {
                'name': name,
                'description': _first_doc_line(func),
                'tags': ['3D'] if _is_3d(func) else ['2D']
            }
            if module_name in _descriptors_tree:
                _descriptors_tree[module_name]['descriptors'].append(descriptor)
            else:
                _descriptors_tree[module_name] = {
                    'version': func.__dict__.get('version', 'unknown'),
                    'name': _module_names.get(module_name, module_name),
                    'description': _first_doc_line(sys.modules[func.__module__]),
                    'descriptors': [descriptor]}

    return _descriptors_tree


def _get_descriptors_list():
    """
    Gets descriptors list.

    :return: Descriptors list: list(tuple(name, func)).
    """
    global _descriptors_list
    if _descriptors_list is None:
        _descriptors_list = Descriptors._descList + [
            ("PMI1", Descriptors3D.PMI1),
            ("PMI2", Descriptors3D.PMI2),
            ("PMI3", Descriptors3D.PMI3),
            ("NPR1", Descriptors3D.NPR1),
            ("NPR2", Descriptors3D.NPR2),
            ("RadiusOfGyration", Descriptors3D.RadiusOfGyration),
            ("InertialShapeFactor", Descriptors3D.InertialShapeFactor),
            ("Eccentricity", Descriptors3D.Eccentricity),
            ("Asphericity", Descriptors3D.Asphericity),
            ("SpherocityIndex", Descriptors3D.SpherocityIndex),
        ]
    return _descriptors_list


def _get_descriptors_funcs(descriptors):
    """
    Gets descriptors functions.

    :param descriptors: Array of descriptors names.
    :return: Dictionary of functions corresponding to descriptors names.
    """
    descriptors_funcs = {}
    descs = _get_descriptors_list()
    for desc in descs:
        if desc[0] in descriptors:
            descriptors_funcs[desc[0]] = desc[1]
    return descriptors_funcs


def get_descriptors(molecules, descriptors):
    """
    Gets descriptors for each input molecule.

    :param molecules: Array of molecules.
    :param descriptors: Array of descriptors or groups names.
    :return: Dictionary of descriptors values for each molecule.
    """
    tree = get_descriptors_tree()
    expanded = []
    for d in descriptors:
        expanded.extend([descriptor['name'] for descriptor in tree[d]['descriptors']] if d in tree else [d])
    descriptors = expanded
    length = len(molecules)
    values = [np_none(length) for _ in descriptors]
    descriptors_funcs = _get_descriptors_funcs(descriptors)
    needs_3d = not set(descriptors).isdisjoint(get_3d_descriptors_names())
    for n in range(0, length):
        try:
            mol = Chem.MolFromMolBlock(molecules[n]) if ("M  END" in molecules[n]) else Chem.MolFromSmiles(molecules[n])
        except:
            mol = None
        if mol is not None and needs_3d:
            try:
                mol = _add_3d_coordinates(mol)
            except:
                mol = None
        if mol is None:
            continue
        for i in range(0, len(descriptors)):
            try:
                value = descriptors_funcs[descriptors[i]](mol)
                if not np.isnan(value) and not np.isinf(value):
                    values[i][n] = value
                    if descriptors[i] not in _descriptor_types:
                        _descriptor_types[descriptors[i]] = "int" if isinstance(value, int) else "double"
            except:
                pass
    return {descriptors[i]: set_type(values[i], _descriptor_types.get(descriptors[i], "double"))
            for i in range(0, len(descriptors))}


def _add_3d_coordinates(mol):
    """
    Add 3D coordinates into molecule.

    :param mol: Molecule.
    :return: Molecule with 3D coordinates.
    """
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    mol = Chem.RemoveHs(mol)
    return mol


def molecules_to_canonical(molecules):
    """
    Converts molecule into canonical SMILES form.

    :param molecules: Array of molecules.
    :return: Array of molecules in canonical SMILES notation.
    """
    length = len(molecules)
    canonical = np_none(length)
    for n in range(0, length):
        try:
            mol = Chem.MolFromMolBlock(molecules[n]) if ("M  END" in molecules[n]) else Chem.MolFromSmiles(molecules[n])
        except:
            mol = None
        canonical[n] = Chem.MolToSmiles(mol) if mol is not None else ''
    return canonical.tolist()


def molecules_to_inchi(molecules):
    """
    Converts molecules into InChI string.

    :param molecules: Array of molecules.
    :return: Array of corresponding InChI string.
    """
    length = len(molecules)
    inchi = np_none(length)
    for n in range(0, length):
        try:
            mol = Chem.MolFromMolBlock(molecules[n]) if ("M  END" in molecules[n]) else Chem.MolFromSmiles(molecules[n])
        except:
            mol = None
        inchi[n] = '' if mol is None else Chem.inchi.MolToInchi(mol)
    return inchi.tolist()


def molecules_to_inchi_key(molecules):
    """
    Converts molecules into InChI key string.

    :param molecules: Array of molecules.
    :return: Array of corresponding InChI key string.
    """
    length = len(molecules)
    inchi_key = np_none(length)
    for n in range(0, length):
        try:
            mol = Chem.MolFromMolBlock(molecules[n]) if ("M  END" in molecules[n]) else Chem.MolFromSmiles(molecules[n])
        except:
            mol = None
        inchi_key[n] = Chem.inchi.MolToInchiKey(mol) if mol != None else ''
    return inchi_key.tolist()


def inchi_to_inchi_key(inchi):
    """
    Converts molecule InChI string into InChI key string.

    :param inchi: Array of molecules in InChI format.
    :return: Array of corresponding InChI key string.
    """
    length = len(inchi)
    inchi_key = np_none(length)
    for n in range(0, length):
        inchi_key[n] = Chem.inchi.InchiToInchiKey(str(inchi[n])) if inchi[n] != '' and inchi[n] is not None else None
    return inchi_key.tolist()


def inchi_to_smiles(inchi):
    """
    Converts molecule InChI string into SMILES.

    :param inchi: Array of molecules in InChI format.
    :return: Array of molecules in SMILES notation.
    """
    length = len(inchi)
    smiles = np_none(length)
    for n in range(0, length):
        try:
            mol = Chem.inchi.MolFromInchi(str(inchi[n]))
        except:
            mol = None
        smiles[n] = Chem.MolToSmiles(mol) if mol != None else ''
    return smiles.tolist()
