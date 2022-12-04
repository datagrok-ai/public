from io import TextIOWrapper

from rdkit import Chem

import orjson

import click

from click_default_group import DefaultGroup
from rdkit.Chem.rdchem import Mol


def molAddCollection(mol: Mol, name: str, title: str = None) -> str:
    """
    Get and postprocess (atom's CFG, title, e.t.c.) molblock
    :param mol:    Mol molecule structure / object
    :param name:   Monomer name to add to molblock title
    :param title:  title to replace in Chem.MolToMolBlock() string output
    :return:       molblock string
    """
    res: str = Chem.MolToMolBlock(mol, forceV3000=True)  # MolToMolFile

    mb_line_list: list[str] = res.split('\n')
    if title:
        mb_line_list[1] = title

    if name and name not in mb_line_list[1]:
        mb_line_list[1] += '|' + name

    end_bond_idx: int = mb_line_list.index('M  V30 END BOND')
    chirality = [atom.GetChiralTag() for atom in mol.GetAtoms()]
    begin_atom_idx = mb_line_list.index('M  V30 BEGIN ATOM')
    end_atom_idx = mb_line_list.index('M  V30 END ATOM')
    for atom_idx in range(1, end_atom_idx - begin_atom_idx):
        line_idx = begin_atom_idx + atom_idx
        atom_ch = chirality[atom_idx - 1]
        if atom_ch != Chem.rdchem.CHI_UNSPECIFIED:
            mb_line_list[line_idx] += " CFG={0}".format(int(atom_ch))

    steabs: list[int] = [i + 1 for (i, ch) in enumerate(chirality) if ch != Chem.rdchem.CHI_UNSPECIFIED]
    if len(steabs) > 0:
        steabs_str: str = "M  V30 MDLV30/STEABS ATOMS=({count} {list})" \
            .format(count=len(steabs), list=' '.join([str(idx) for idx in steabs]))

        mb_line_list = mb_line_list[:(end_bond_idx + 1)] + \
                       ["M  V30 BEGIN COLLECTION", steabs_str, "M  V30 END COLLECTION"] + \
                       mb_line_list[(end_bond_idx + 1):]

    return '\n'.join(mb_line_list)


def molfile2molfile(src_mol: str, name: str) -> str:
    mol: Mol = Chem.MolFromMolBlock(src_mol)
    src_mf_lines = src_mol.split('\n')
    title = src_mf_lines[1]
    return molAddCollection(mol, name, title=title)


def smiles2molfile(smiles: str, name: str) -> str:
    mol: Mol = Chem.MolFromSmiles(smiles)
    return molAddCollection(mol, name)


CodesType = dict[str, dict[str, list[str]]]


class Monomer:
    def __init__(self,
                 symbol: str, name: str, molfile: str, smiles: str,
                 codes: CodesType):
        self.monomerType = 'Backbone'
        self.smiles = smiles
        self.name = name
        self.author = 'SequenceTranslator'
        self.molfile = molfile2molfile(molfile, name) if molfile else smiles2molfile(smiles, name)
        self.naturalAnalog = ''
        self.rgroups = [
            {
                "capGroupSmiles": "O[*:1]",
                "alternateId": "R1-OH",
                "capGroupName": "OH",
                "label": "R1"
            },
            {
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
                    monomers_res[monomer_name] = Monomer(symbol, name, None, smiles, {})
                codes = monomers_res[monomer_name].codes
                if codes_src not in codes:
                    codes[codes_src] = {}
                if codes_type not in codes[codes_src]:
                    codes[codes_src][codes_type] = [];
                codes[codes_src][codes_type].append(codes_code)
    return monomers_res


@click.group(cls=DefaultGroup, default='main')
def cli():
    pass


@cli.command()
@click.pass_context
@click.option('--initial', 'initial_f',
              help='Initial monomers source file.',
              type=click.File('r', 'utf-8'))
@click.option('--lib', 'lib_f',
              help='Output library (HELM format) file.',
              type=click.File('wb', 'utf-8'))
@click.option('--add', 'add_f_list', multiple=True,
              help='Additional libraries to build.',
              type=click.File('r', 'utf-8'))
def main(ctx, initial_f: TextIOWrapper, lib_f: TextIOWrapper, add_f_list: list[TextIOWrapper]):
    initial_json_str = initial_f.read()

    initial_json = orjson.loads(initial_json_str)

    monomers: dict[str, Monomer] = codes2monomers(initial_json)

    for add_f in add_f_list:
        add_json_str = add_f.read()
        add_json = orjson.loads(add_json_str)
        for add_m in add_json:
            m = Monomer.from_json(add_m)
            monomers[m.name] = m

    add_json = [m.to_json() for m in monomers.values()]

    lib_json_txt = orjson.dumps(add_json, option=orjson.OPT_INDENT_2)
    lib_f.write(lib_json_txt)
    k = 11


if __name__ == '__main__':
    cli()
