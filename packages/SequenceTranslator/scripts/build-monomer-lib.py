# pylint: disable=no-member
import os.path
import sys
from io import TextIOWrapper
from typing import Optional

from rdkit import Chem
from rdkit.Chem.rdchem import Mol

import orjson

import click

CodesType = dict[str, dict[str, list[str]]]

BEGIN_ATOM_LINE = 'M  V30 BEGIN ATOM'
END_ATOM_LINE = 'M  V30 END ATOM'
BEGIN_BOND_LINE = 'M  V30 BEGIN BOND'
END_BOND_LINE = 'M  V30 END BOND'
BEGIN_COLLECTION_LINE = 'M  V30 BEGIN COLLECTION'
END_COLLECTION_LINE = 'M  V30 END COLLECTION'
COLLECTION_STEABS_LINE = 'M  V30 MDLV30/STEABS'
IDX_OF_FIRST_VALUE = 7
NUM_OF_BOND_POSITIONAL_ARGS = 4
CFG = "CFG="


def mol_add_collection(mol: Mol,
                       name: str,
                       title: Optional[str] = None,
                       src_mol: Optional[str] = None) -> str:
    """
    Get and postprocess (atom's CFG, title, e.t.c.) molblock
    :param mol:      Mol molecule structure / object
    :param name:     Monomer name to add to molblock title string
    :param title:    Title to replace in Chem.MolToMolBlock() string output
    :param src_mol:  Source molblock data, to restore optional CFG
    :return:         molblock string
    """
    res: str = Chem.MolToMolBlock(mol, forceV3000=True)  # MolToMolFile

    mb_line_list: list[str] = res.split('\n')
    if title:
        mb_line_list[1] = title

    if name and name not in mb_line_list[1]:
        mb_line_list[1] += '|' + name

    chirality = [atom.GetChiralTag() for atom in mol.GetAtoms()]

    # preserve chirality for bonds from src_mol
    tgt_mol_file_map = MolFileMap.parse(res)
    steabs = []
    if src_mol:
        src_mol_file_map = MolFileMap.parse(src_mol)
        if len(tgt_mol_file_map.mol_file.atom_list) != len(src_mol_file_map.mol_file.atom_list):
            raise ValueError(f"Atoms count of src and tgt differs for monomer '{name}'.")
        for (src_bond_idx0, (bond_key, src_bond)) in enumerate(src_mol_file_map.mol_file.bonds.items()):
            if src_bond.cfg:
                if bond_key not in tgt_mol_file_map.mol_file.bonds:
                    raise KeyError(f"Bond key '{bond_key}' not found in tgt bonds.")
                tgt_bond: MolFileBond = tgt_mol_file_map.mol_file.bonds[bond_key]
                tgt_bond_cfg_str: str = ' '.join(tgt_bond.cfg)
                src_bond_cfg_str: str = ' '.join(src_bond.cfg)
                if tgt_bond_cfg_str != src_bond_cfg_str:
                    mb_line_list[tgt_mol_file_map.begin_bond_idx + tgt_bond.bond_idx] += f" {src_bond_cfg_str}"

        for (tgt_atom_idx0, tgt_atom) in enumerate(tgt_mol_file_map.mol_file.atom_list):
            src_atom = src_mol_file_map.mol_file.atom_list[tgt_atom_idx0]
            atom_ch = chirality[tgt_atom_idx0]
            if src_atom.cfg:
                mb_line_list[tgt_mol_file_map.begin_atom_idx + tgt_atom_idx0 + 1] += " {0}".format(
                    ' '.join(src_atom.cfg))
                steabs.append(tgt_atom_idx0 + 1)
            elif atom_ch != Chem.rdchem.CHI_UNSPECIFIED:
                mb_line_list[tgt_mol_file_map.begin_atom_idx + tgt_atom_idx0 + 1] += " CFG={0}".format(int(atom_ch))
                steabs.append(tgt_atom_idx0 + 1)
            elif src_atom.atom_idx in src_mol_file_map.mol_file.collection_steabs:
                raise KeyError(f"Source STEABS atom '{src_atom}' not accounted")
            elif tgt_atom.atom_idx in tgt_mol_file_map.mol_file.collection_steabs:
                raise KeyError(f"Target STEABS atom '{tgt_atom}' not accounted")
    else:
        # from Smiles (without src_mol)
        for (tgt_atom_idx0, tgt_atom) in enumerate(tgt_mol_file_map.mol_file.atom_list):
            line_idx = tgt_mol_file_map.begin_atom_idx + tgt_atom_idx0 + 1
            atom_ch = chirality[tgt_atom_idx0]
            if atom_ch != Chem.rdchem.CHI_UNSPECIFIED:
                mb_line_list[line_idx] += " CFG={0}".format(int(atom_ch))
                steabs.append(tgt_atom_idx0 + 1)

    if len(steabs) > 0:
        steabs_str: str = COLLECTION_STEABS_LINE + " ATOMS=({count} {list})".format(
            count=len(steabs),
            list=' '.join([str(idx) for idx in steabs]))
        if tgt_mol_file_map.collection_steabs_idx:
            mb_line_list[tgt_mol_file_map.collection_steabs_idx] = steabs_str
        elif tgt_mol_file_map.begin_collection_idx is not None:
            tgt_collection_steabs_idx = tgt_mol_file_map.begin_collection_idx + 1
            mb_line_list = mb_line_list[:tgt_collection_steabs_idx] + \
                           [steabs_str] + \
                           mb_line_list[tgt_collection_steabs_idx:]
        else:
            tgt_collection_idx = tgt_mol_file_map.end_bond_idx + 1
            mb_line_list = mb_line_list[:tgt_collection_idx] + \
                           [BEGIN_COLLECTION_LINE, steabs_str, END_COLLECTION_LINE] + \
                           mb_line_list[tgt_collection_idx:]

    return '\n'.join(mb_line_list)


def prepare_molblock(src_molblock: str, name: str) -> str:
    """Loads mol from src_mol str. Fixed title, adds chirality to atoms and preserves chirality for bonds."""
    # Using sanitize=False leads to unwanted moving stereo (invalid?) CFGs to other bonds
    mol: Mol = Chem.MolFromMolBlock(src_molblock, removeHs=False)
    src_molblock_lines = src_molblock.split('\n')
    title = src_molblock_lines[1]
    return mol_add_collection(mol, name, title=title, src_mol=src_molblock)


def smiles2molfile(smiles: str, name: str) -> str:
    mol: Mol = Chem.MolFromSmiles(smiles)
    return mol_add_collection(mol, name)


class Monomer:
    def __init__(self, symbol: str, name: str, molfile: str, smiles: str,
                 codes: CodesType):
        self.monomerType = 'Backbone'
        self.smiles = smiles
        self.name = name
        self.author = 'SequenceTranslator'
        self.molfile = prepare_molblock(
            molfile, name) if molfile else smiles2molfile(smiles, name)
        self.naturalAnalog = ''
        self.rgroups = [{
            "capGroupSmiles": "O[*:1]",
            "alternateId": "R1-OH",
            "capGroupName": "OH",
            "label": "R1"
        }, {
            "capGroupSmiles": "O[*:2]",
            "alternateId": "R2-OH",
            "capGroupName": "OH",
            "label": "R2"
        }]
        self.createDate = None
        self.id = 0
        self.polymerType = 'RNA'
        self.symbol = symbol
        self.codes: CodesType = codes

    @staticmethod
    def from_json(src_json: {}):
        obj = Monomer(src_json['symbol'], src_json['name'],
                      src_json['molfile'], src_json['smiles'],
                      src_json['codes'])
        return obj

    def to_json(self):
        return {
            'monomerType': self.monomerType,
            'smiles': self.smiles,
            'name': self.name,
            'author': self.author,
            'molfile': self.molfile,
            'naturalAnalog': self.naturalAnalog,
            'rgroups': self.rgroups,
            'createDate': self.createDate,
            'id': self.id,
            'polymerType': self.polymerType,
            'symbol': self.symbol,
            'codes': self.codes,
        }


class MolFileAtom:
    """
    Wrapper for data extracted from molfile atom line
    """
    def __init__(self, v3k_atom_line: str):
        self._atom_line = v3k_atom_line
        self._atom_line_splitted: [] = self.\
            _atom_line[IDX_OF_FIRST_VALUE:].split(' ')
        self._atom_idx = int(self._atom_line_splitted[0].strip())
        # we cannot use positional argument for cfg for it is a kwarg
        cfg_item = list(filter(
            lambda x: x.startswith(CFG), self._atom_line_splitted
            ))
        self._cfg = cfg_item

    @property
    def atom_line_str(self) -> str:
        return self._atom_line

    @property
    def atom_idx(self):
        return self._atom_idx

    @property
    def atom_line_splitted(self) -> list[str]:
        return self.atom_line_splitted

    @property
    def cfg(self) -> list[str]:
        return self._cfg

    @property
    def cfg_int(self) -> int:
        return self._cfg

    def __str__(self):
        return self._atom_line

    def __repr__(self):
        return str(self)


class MolFileBond:
    """
    Wrapper for data extracted from molfile bond line
    """
    def __init__(self, v3k_bond_line: str):
        self._bond_line = v3k_bond_line
        self._bond_line_splitted: [] = self.\
            _bond_line[IDX_OF_FIRST_VALUE:].split(' ')
        self._bond_idx = int(self._bond_line_splitted[0].strip())
        self._key = self._bond_line_splitted[0:NUM_OF_BOND_POSITIONAL_ARGS]
        cfg_item = list(filter(
            lambda x: x.startswith(CFG), self._bond_line_splitted
            ))
        self._cfg = cfg_item

    @property
    def bond_line(self) -> str:
        return self._bond_line

    @property
    def bond_idx(self):
        return self._bond_idx

    @property
    def bond_line_splitted(self) -> list[str]:
        return self._bond_line_splitted

    @property
    def key(self):
        return self._key

    @property
    def cfg(self) -> list[str]:
        return self._cfg

    def __str__(self):
        return self._bond_line

    def __repr__(self):
        return str(self)


class MolFileV3K:
    """
    Wrapper for data extracted from molfile
    """
    def __init__(
            self, title: str, atom_list: list[MolFileAtom],
            bond_list: list[MolFileBond],
            collection_steabs: list[int] = None
            ):
        self._title = title
        self._atom_list = atom_list
        self._bond_list = bond_list
        self.collection_steabs = [] if collection_steabs is None \
            else collection_steabs

        self._bonds: dict = {}
        for bond in self._bond_list:
            # list is unhashable type, but tuple is
            bond_key = tuple((int(v) for v in bond.key))
            self._bonds[bond_key] = bond

    @property
    def atom_list(self):
        return self._atom_list

    @property
    def bond_list(self):
        return self._bond_list

    @property
    def bonds(self) -> dict:
        return self._bonds


class MolFileMap:
    def __init__(self, src: str, mol_file_obj: MolFileV3K,
                 atom_block_idx_boundaries: tuple[int, int],
                 bond_block_idx_boundaries: tuple[int, int],
                 collection_idx_boundaries: tuple[int, int] = None,
                 collection_steabs_idx: int = None):
        self._src = src
        self._mol_file = mol_file_obj
        self.begin_atom_idx = atom_block_idx_boundaries[0]
        self.end_atom_idx = atom_block_idx_boundaries[1]
        self.begin_bond_idx = bond_block_idx_boundaries[0]
        self.end_bond_idx = bond_block_idx_boundaries[1]
        self.begin_collection_idx = None if collection_idx_boundaries is None \
            else collection_idx_boundaries[0]
        self.end_collection_idx = None if collection_idx_boundaries is None \
            else collection_idx_boundaries[1]
        self.collection_steabs_idx = collection_steabs_idx

    @property
    def src(self):
        return self._src

    @property
    def mol_file(self):
        return self._mol_file

    @staticmethod
    def parse(molblock_src: str):
        molblock_line_list: list[str] = \
            [line.rstrip() for line in molblock_src.split('\n')]
        title: str = molblock_line_list[1]

        def get_idx_boundaries(begin_str: str, end_str: str):
            return tuple([
                molblock_line_list.index(begin_str),
                molblock_line_list.index(end_str)
                ])

        def get_wrapper_list(
                begin_idx: int, end_idx: int, wrapper_constructor
                ):
            """
            For the list of atom/bond wrapper objects
            """
            item_count = end_idx - begin_idx - 1  # for atoms or bonds
            wrapper_list = [None] * item_count
            for item_idx in range(1, item_count + 1):
                line_idx = begin_idx + item_idx
                line = molblock_line_list[line_idx]
                item = wrapper_constructor(line)
                wrapper_list[item_idx - 1] = item
            return wrapper_list

        atom_block_idx_boundaries = get_idx_boundaries(
                BEGIN_ATOM_LINE, END_ATOM_LINE)
        bond_block_idx_boundaries = get_idx_boundaries(
                BEGIN_BOND_LINE, END_BOND_LINE)
        atom_list = get_wrapper_list(
                atom_block_idx_boundaries[0],
                atom_block_idx_boundaries[1], MolFileAtom)
        bond_list = get_wrapper_list(
                bond_block_idx_boundaries[0],
                bond_block_idx_boundaries[1], MolFileBond)

        collection_idx_boundaries = None
        collection_steabs_idx = None
        collection_steabs: list[int] = []
        if BEGIN_COLLECTION_LINE in molblock_line_list and END_COLLECTION_LINE in molblock_line_list:
            collection_idx_boundaries = get_idx_boundaries(
                    BEGIN_COLLECTION_LINE, END_COLLECTION_LINE)
            collection_count: int = collection_idx_boundaries[1] - \
                collection_idx_boundaries[0] - 1
            for collection_idx in range(1, collection_count + 1):
                line_idx = collection_idx_boundaries[0] + collection_idx
                collection_line = molblock_line_list[line_idx]
                if collection_line.startswith(COLLECTION_STEABS_LINE):
                    steabs_str = collection_line[len(COLLECTION_STEABS_LINE + " ATOMS=("):-1]
                    collection_steabs = [int(atom_num_str.strip()) for atom_num_str in steabs_str.split(' ')[1:]]
                    collection_steabs_idx = line_idx
                else:
                    raise ValueError(f"Unexpected collection line '{collection_line}'.")

        mol_file = MolFileV3K(title, atom_list, bond_list, collection_steabs)
        return MolFileMap(
                molblock_src, mol_file,
                atom_block_idx_boundaries, bond_block_idx_boundaries,
                collection_idx_boundaries, collection_steabs_idx)


def codes2monomers(codes_json: {}) -> dict[str, Monomer]:
    resulting_monomer_dict: dict[str, Monomer] = {}
    for (codes_src, src_dict) in codes_json.items():
        for (codes_type, monomers_dict) in src_dict.items():
            for (codes_code, monomer_json) in monomers_dict.items():
                monomer_name = monomer_json['name']
                if monomer_name not in resulting_monomer_dict:
                    symbol = monomer_json['name']
                    name = monomer_json['name']
                    smiles = monomer_json['SMILES']
                    resulting_monomer_dict[monomer_name] = Monomer(
                        symbol, name, None, smiles, {})
                codes = resulting_monomer_dict[monomer_name].codes
                if codes_src not in codes:
                    codes[codes_src] = {}
                if codes_type not in codes[codes_src]:
                    codes[codes_src][codes_type] = []
                codes[codes_src][codes_type].append(codes_code)
    return resulting_monomer_dict


@click.command()
@click.option('--initial',
              'initial_f',
              help='Initial monomers source file.',
              type=click.File('r', 'utf-8'))
@click.option('--lib',
              'lib_f',
              help='Output library (HELM format) file.',
              type=click.File('wb', 'utf-8'))
@click.option('--add-list',
              'monomer_list_file_list',
              multiple=True,
              help='Monomers to be added to the library',
              type=click.File('r', 'utf-8'))
def main(initial_f: TextIOWrapper, lib_f: TextIOWrapper,
         monomer_list_file_list: list[TextIOWrapper]):
    monomers: dict[str, Monomer] = {}
    if initial_f:
        initial_json_str = initial_f.read()
        initial_json = orjson.loads(initial_json_str)
        monomers.update(codes2monomers(initial_json))

    monomer_filename_list = []
    for monomer_list_file in monomer_list_file_list:
        for monomer_fn in [fn for fn in monomer_list_file.read().split('\n') if fn]:
            monomer_filename_list.append(monomer_fn);

    print(monomer_filename_list)

    for add_json_fn in monomer_filename_list:
        # trying to load mol data if file with .mol extension exists
        with open(add_json_fn, 'r') as add_json_f:
            add_json_str = add_json_f.read()
            add_json = orjson.loads(add_json_str)
            for add_m in add_json:
                name = add_m['name']
                try:
                    if not add_m['molfile']:
                        monomer_mol_fn = os.path.join(os.path.dirname(add_json_fn), add_m['name'] + '.mol')
                        if not os.path.isfile(monomer_mol_fn):
                            raise FileNotFoundError(monomer_mol_fn)
                        else:
                            with open(monomer_mol_fn, 'r') as monomer_mol_f:
                                monomer_mol_lines = [line.rstrip() for line in monomer_mol_f.readlines()]
                                monomer_mol_txt = '\n'.join(monomer_mol_lines)
                                add_m['molfile'] = monomer_mol_txt

                    m = Monomer.from_json(add_m)
                    monomers[m.name] = m
                except Exception as ex:
                    sys.stderr.write(f"Invalid monomer '{add_m['name']}' error:\n{str(ex)}")

    add_json = [m.to_json() for m in monomers.values()]

    lib_json_txt = orjson.dumps(add_json, option=orjson.OPT_INDENT_2)
    lib_f.write(lib_json_txt)


if __name__ == '__main__':
    main()
