# pylint: disable=no-member
import os.path
import sys
from io import TextIOWrapper
from typing import Optional

from rdkit import Chem
from rdkit.Chem.rdchem import Mol

import orjson

import click

from click_default_group import DefaultGroup

BEGIN_ATOM_LINE = 'M  V30 BEGIN ATOM'
END_ATOM_LINE = 'M  V30 END ATOM'
BEGIN_BOND_LINE = 'M  V30 BEGIN BOND'
END_BOND_LINE = 'M  V30 END BOND'
BEGIN_COLLECTION_LINE = 'M  V30 BEGIN COLLECTION'
END_COLLECTION_LINE = 'M  V30 END COLLECTION'
COLLECTION_STEABS_LINE = 'M  V30 MDLV30/STEABS'


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
                    mb_line_list[tgt_mol_file_map.begin_bond_idx + tgt_bond.number] += f" {src_bond_cfg_str}"

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
            elif src_atom.number in src_mol_file_map.mol_file.collection_steabs:
                raise KeyError(f"Source STEABS atom '{src_atom}' not accounted")
            elif tgt_atom.number in tgt_mol_file_map.mol_file.collection_steabs:
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


def molfile2molfile(src_mol: str, name: str) -> str:
    """Loads mol from src_mol str. Fixed title, adds chirality to atoms and preserves chirality for bonds."""
    # Using sanitize=False leads to unwanted moving stereo (invalid?) CFGs to other bonds
    mol: Mol = Chem.MolFromMolBlock(src_mol, removeHs=False)
    src_mf_lines = src_mol.split('\n')
    title = src_mf_lines[1]
    return mol_add_collection(mol, name, title=title, src_mol=src_mol)


def smiles2molfile(smiles: str, name: str) -> str:
    mol: Mol = Chem.MolFromSmiles(smiles)
    return mol_add_collection(mol, name)


CodesType = dict[str, dict[str, list[str]]]


class Monomer:

    def __init__(self, symbol: str, name: str, molfile: str, smiles: str,
                 codes: CodesType):
        self.monomerType = 'Backbone'
        self.smiles = smiles
        self.name = name
        self.author = 'SequenceTranslator'
        self.molfile = molfile2molfile(
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
    def __init__(self, src: str):
        self._src = src
        self._parts: [] = self._src[7:].split(' ')
        self._number = int(self._parts[0].strip())
        self._cfg = self._parts[6:]

    @property
    def src(self) -> str:
        return self._src

    @property
    def number(self):
        return self._number

    @property
    def parts(self) -> list[str]:
        return self.parts

    @property
    def cfg(self) -> list[str]:
        return self._cfg

    def __str__(self):
        return self._src

    def __repr__(self):
        return str(self)


class MolFileBond:
    def __init__(self, src: str):
        self._src = src
        self._parts: [] = self._src[7:].split(' ')
        self._number = int(self._parts[0].strip())
        self._key = self._parts[0:4]
        self._cfg = self._parts[4:]

    @property
    def src(self) -> str:
        return self._src

    @property
    def number(self):
        return self._number

    @property
    def parts(self) -> list[str]:
        return self.parts

    @property
    def key(self):
        return self._key

    @property
    def cfg(self) -> list[str]:
        return self._cfg

    def __str__(self):
        return self._src

    def __repr__(self):
        return str(self)


class MolFile:
    def __init__(self, title: str, atom_list: list[MolFileAtom], bond_list: list[MolFileBond],
                 collection_steabs: list[int] = []):
        self._title = title
        self._atom_list = atom_list
        self._bond_list = bond_list
        self.collection_steabs = collection_steabs

        self._bonds: dict = {}
        for bond in self._bond_list:
            bond_line = bond.src
            bond_parts = bond_line[7:].split(' ')
            # list is unhashable type, but tuple s
            bond_key = tuple((int(v) for v in bond_parts[0:4]))
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
    def __init__(self, src: str, mol_file: MolFile,
                 begin_atom_idx: int, end_atom_idx: int, begin_bond_idx: int, end_bond_idx: int,
                 begin_collection_idx: int = None, end_collection_idx: int = None,
                 collection_steabs_idx: int = None):
        self._src = src
        self._mol_file = mol_file
        self.begin_atom_idx = begin_atom_idx
        self.end_atom_idx = end_atom_idx
        self.begin_bond_idx = begin_bond_idx
        self.end_bond_idx = end_bond_idx
        self.begin_collection_idx = begin_collection_idx
        self.end_collection_idx = end_collection_idx
        self.collection_steabs_idx = collection_steabs_idx

    @property
    def src(self):
        return self._src

    @property
    def mol_file(self):
        return self._mol_file

    @staticmethod
    def parse(mol_src: str):
        mb_line_list: list[str] = [line.rstrip() for line in mol_src.split('\n')]
        title: str = mb_line_list[1]

        begin_atom_idx: int = mb_line_list.index(BEGIN_ATOM_LINE)
        end_atom_idx: int = mb_line_list.index(END_ATOM_LINE)
        atom_count = end_atom_idx - begin_atom_idx - 1
        atom_list: list[Optional[MolFileAtom]] = [None] * atom_count
        for atom_idx in range(1, atom_count + 1):
            line_idx = begin_atom_idx + atom_idx
            atom_line = mb_line_list[line_idx]
            atom = MolFileAtom(atom_line)
            atom_list[atom_idx - 1] = atom

        begin_bond_idx: int = mb_line_list.index(BEGIN_BOND_LINE)
        end_bond_idx: int = mb_line_list.index(END_BOND_LINE)
        bond_count: int = end_bond_idx - begin_bond_idx - 1
        bond_list: list[Optional[MolFileBond]] = [None] * bond_count
        for bond_idx in range(1, bond_count + 1):
            line_idx = begin_bond_idx + bond_idx
            bond_line = mb_line_list[line_idx]
            bond = MolFileBond(bond_line)
            bond_list[bond_idx - 1] = bond

        begin_collection_idx = None
        end_collection_idx = None
        collection_steabs_idx = None
        collection_steabs: list[int] = []
        if BEGIN_COLLECTION_LINE in mb_line_list and END_COLLECTION_LINE in mb_line_list:
            begin_collection_idx = mb_line_list.index(BEGIN_COLLECTION_LINE)
            end_collection_idx = mb_line_list.index(END_COLLECTION_LINE)
            collection_count: int = end_collection_idx - begin_collection_idx - 1
            for collection_idx in range(1, collection_count + 1):
                line_idx = begin_collection_idx + collection_idx
                collection_line = mb_line_list[line_idx]
                if collection_line.startswith(COLLECTION_STEABS_LINE):
                    steabs_str = collection_line[len(COLLECTION_STEABS_LINE + " ATOMS=("):-1]
                    collection_steabs = [int(atom_num_str.strip()) for atom_num_str in steabs_str.split(' ')[1:]]
                    collection_steabs_idx = line_idx
                else:
                    raise ValueError(f"Unexpected collection line '{collection_line}'.")

        mol_file = MolFile(title, atom_list, bond_list, collection_steabs)
        return MolFileMap(mol_src, mol_file,
                          begin_atom_idx, end_atom_idx, begin_bond_idx, end_bond_idx,
                          begin_collection_idx, end_collection_idx,
                          collection_steabs_idx)


def codes2monomers(codes_json: {}) -> dict[str, Monomer]:
    monomers_res: dict[str, Monomer] = {}
    for (codes_src, src_dict) in codes_json.items():
        for (codes_type, monomers_dict) in src_dict.items():
            for (codes_code, monomer_json) in monomers_dict.items():
                monomer_name = monomer_json['name']
                if monomer_name not in monomers_res:
                    symbol = monomer_json['name']
                    name = monomer_json['name']
                    smiles = monomer_json['SMILES']
                    monomers_res[monomer_name] = Monomer(
                        symbol, name, None, smiles, {})
                codes = monomers_res[monomer_name].codes
                if codes_src not in codes:
                    codes[codes_src] = {}
                if codes_type not in codes[codes_src]:
                    codes[codes_src][codes_type] = []
                codes[codes_src][codes_type].append(codes_code)
    return monomers_res


@click.group(cls=DefaultGroup, default='main')
def cli():
    pass


@cli.command()
@click.pass_context
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
              help='Additional libraries to build.',
              type=click.File('r', 'utf-8'))
def main(ctx, initial_f: TextIOWrapper, lib_f: TextIOWrapper,
         monomer_list_file_list: list[TextIOWrapper]):
    monomers: dict[str, Monomer] = {}
    if initial_f:
        initial_json_str = initial_f.read()
        initial_json = orjson.loads(initial_json_str)
        monomers.update(codes2monomers(initial_json))

    monomer_fn_list = []
    for monomer_list_file in monomer_list_file_list:
        for monomer_fn in [fn for fn in monomer_list_file.read().split('\n') if fn]:
            monomer_fn_list.append(monomer_fn);

    print(monomer_fn_list)

    for add_json_fn in monomer_fn_list:
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
    cli()
