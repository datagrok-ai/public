from io import TextIOWrapper

from rdkit.Chem import AllChem
from rdkit import Chem

import orjson
import json

import click

from click_default_group import DefaultGroup
from rdkit.Chem.rdchem import Mol


def smiles2molfile(smiles: str) -> str:
    mol: Mol = Chem.MolFromSmiles(smiles)
    res: str = Chem.MolToMolBlock(mol, forceV3000=True)  # MolToMolFile
    return res


def molV2000toMolV3000(molV2K: str) -> str:
    mol: str = Chem.MolFromMolBlock(molV2K)
    res: str = Chem.MolToMolBlock(mol, forceV3000=True)
    return res.replace('Pol', 'O  ')


CodesType = dict[str, dict[str, list[str]]]


class Monomer:
    def __init__(self,
                 symbol: str, name: str, smiles: str,
                 codes: CodesType):
        self.monomerType = 'Backbone'
        self.smiles = smiles
        self.name = name
        self.author = 'SequenceTranslator'
        self.molfile = smiles2molfile(smiles)
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
        obj = Monomer(src_json['symbol'], src_json['name'], src_json['smiles'], src_json['codes'])
        obj.molfile = src_json['molfile']
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
                    monomers_res[monomer_name] = Monomer(symbol, name, smiles, {})
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
