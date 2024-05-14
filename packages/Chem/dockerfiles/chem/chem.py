from rdkit.Chem import Descriptors
from rdkit.Chem import Descriptors3D
from rdkit.Chem import AllChem
import rdkit.Chem as Chem
import numpy as np

from chem_utils import get_3d_descriptors_names, get_descriptor_description, get_descriptor_tags, get_descriptor_type, \
    get_module_description, get_module_name
from utils import np_none, set_type

_descriptors_tree = None


def get_descriptors_tree():
    """
    Gets all available Molecule descriptors in tree view (module -> descriptor).

    :return: List of available Molecule descriptors.
    """
    global _descriptors_tree
    if _descriptors_tree is None:
        _descriptors_tree = {}
        for desc in _get_descriptors_list():
            func = desc[1]
            name = desc[0]
            module_name = func.__module__.replace('rdkit.Chem.', '')
            version = func.__dict__['version'] if ('version' in func.__dict__) else 'unknown'
            descriptor = {
                'name': name,
                'description': get_descriptor_description(name),
                'type': get_descriptor_tags(name),
                'tags': get_descriptor_tags(name)
            }
            if module_name in _descriptors_tree:
                _descriptors_tree[module_name]['descriptors'].append(descriptor)
            else:
                _descriptors_tree[module_name] = {
                    'version': version,
                    'name': get_module_name(module_name),
                    'description': get_module_description(module_name),
                    'descriptors': [descriptor]}

    return _descriptors_tree


def _get_descriptors_list():
    """
    Gets descriptors list.

    :return: Descriptors list: list(tuple(name, func)).
    """
    descs = Descriptors._descList[:]
    descs_extra = [
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
    for de in descs_extra:
        descs.append(de)
    return descs


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
    _descriptors = list(descriptors)
    tree = get_descriptors_tree()
    for group in tree:
        if group in _descriptors:
            descriptors.remove(group)
            for descriptor in tree[group]['descriptors']:
                descriptors.append(descriptor['name'])
    length = len(molecules)
    values = []
    for _ in descriptors:
        values.append(np_none(length))
    descriptors_3d = get_3d_descriptors_names()
    descriptors_funcs = _get_descriptors_funcs(descriptors)
    for n in range(0, length):
        try:
            mol = Chem.MolFromMolBlock(molecules[n]) if ("M  END" in molecules[n]) else Chem.MolFromSmiles(molecules[n])
        except:
            mol = ''
        try:
            for d in descriptors:
                if d in descriptors_3d and mol != '':
                    try:
                        mol = _add_3d_coordinates(mol)
                    except:
                        mol = ''
                    break
            for d in range(0, len(descriptors)):
                value = descriptors_funcs[descriptors[d]](mol)
                if not np.isnan(value) and not np.isinf(value):
                    values[d][n] = value if get_descriptor_type(descriptors[d]) != "string" else str(value)
        except:
            values[d][n] = None
    result = {}
    for d in range(0, len(descriptors)):
        result[descriptors[d]] = set_type(values[d], get_descriptor_type(descriptors[d]))

    return result


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
