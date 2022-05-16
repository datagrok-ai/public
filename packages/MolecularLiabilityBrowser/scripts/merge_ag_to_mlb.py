from types import SimpleNamespace

import pandas as pd
import orjson as json

from tqdm import tqdm
from numpy import nan, argmax
from string import ascii_letters
import re

import click
from click_default_group import DefaultGroup


@click.group(cls=DefaultGroup, default='main')
def cli():
    pass


@cli.command()
@click.option('--mlb', 'mlb_f',
              help="CSV file with original mlb table with 'v id' column",
              type=click.File('r'))
@click.option('--ab2ag', 'ab2ag_f',
              help="CSV file with antibody to antigen data file with columns 'v_id', 'antigen_list'",
              type=click.File('r'))
@click.option('--ag2g', 'ag2g_f',
              help="CSV file antigen to gene data file with columns 'antigen', 'antigen_ncbi_id', 'antigen_gene_symbol'",
              type=click.File('r'))
@click.option('--out', 'out_f',
              help="CSV file output with joined data",
              type=click.File('w'))
@click.pass_context
def main(ctx, mlb_f, ab2ag_f, ag2g_f, out_f):
    mlb_df: pd.DataFrame = pd.read_csv(mlb_f, dtype=str)

    if 'antigen_list' in mlb_df or 'antigen_ncbi_id' in mlb_df or 'antigen_gene_symbol' in mlb_df:
        raise ValueError("Source MLB file already contains antigen data")

    ab2ag_df = pd.read_csv(ab2ag_f)
    ag2g_df = pd.read_csv(ag2g_f)

    ab2ag_df.antigen_list = [
        re.sub(r'[{}]', '', v) if v != '[Not Found]' else None
        for v in ab2ag_df.antigen_list]

    out_df = pd.merge(mlb_df, ab2ag_df.rename({'v_id': 'v id'}, axis=1), how='left', left_on='v id', right_on='v id')

    ag_dict = dict((
        (row.antigen, SimpleNamespace(ncbi_id=str(row.antigen_ncbi_id), gene_symbol=row.antigen_gene_symbol))
        for (_, row) in ag2g_df.iterrows()))

    # ab2ag data have a strange structure with lists in column 'antigen_list'
    def adds(row) -> (str, str):
        antigen_list_str = row.antigen_list
        antigen_list = antigen_list_str.split(',') if antigen_list_str is not None else []

        ncbi_id_set = set()
        gene_symbol_set = set()
        for antigen in antigen_list:
            ncbi_id_set.add(ag_dict[antigen].ncbi_id)
            gene_symbol_set.add(ag_dict[antigen].gene_symbol)
        return (
            ','.join(ncbi_id_set) if len(ncbi_id_set) > 0 else None,
            ','.join(gene_symbol_set) if len(gene_symbol_set) > 0 else None)

    out_df[['antigen_ncbi_id', 'antigen_gene_symbol']] = out_df.apply(adds, axis=1, result_type='expand')

    out_df.to_csv(out_f, index=False, quoting=False, line_terminator='\n')

    k = 11


if __name__ == '__main__':
    cli()
