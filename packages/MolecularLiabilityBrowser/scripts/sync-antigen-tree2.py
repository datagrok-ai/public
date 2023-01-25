#!/usr/bin/env python

"""
Lists directories for antigen igphyml results.
Imports data from .tab file to mlb.tree table.
"""
import re
import sys
import traceback

antigen_results_re = re.compile('AB_AG_DATA_([A-Z0-9]+)')

import os.path
from datetime import datetime
from types import SimpleNamespace

import numpy as np
import pandas as pd

import click
from click_default_group import DefaultGroup

import sqlalchemy as sa
import sqlalchemy.dialects.postgresql as sapg
import sqlalchemy.engine as sae
import sqlalchemy.orm as sao

meta_public = sa.MetaData(schema='public')
meta_v2 = sa.MetaData(schema='db_v2')
meta_mlb = sa.MetaData(schema='mlb')


def build_schema():
    antigen = sa.Table(
        'antigen', meta_mlb,
        sa.Column('id', sa.Integer, primary_key=True, autoincrement=True),
        sa.Column('antigen', sa.String(), nullable=False, unique=True),
        sa.Column('antigen_ncbi_id', sa.String(), nullable=True),
        sa.Column('antigen_gene_symbol', sa.String(), nullable=True)
    )

    tree2 = sa.Table(
        'tree2', meta_mlb,
        sa.Column('tree_id', sa.Integer, primary_key=True, autoincrement=True),
        sa.Column('antigen_id', sa.Integer, sa.ForeignKey(antigen.columns['id']), nullable=False),
        sa.Column('antigen', sa.String(), sa.ForeignKey(antigen.columns['antigen']), nullable=False),
        sa.Column('CLONE', sa.String(), nullable=False),  # also for 'REPERTOIR'
        sa.Column('NSEQ', sa.Integer, nullable=False),
        sa.Column('NSITE', sapg.DOUBLE_PRECISION, nullable=False),
        sa.Column('TREE_LENGTH', sapg.DOUBLE_PRECISION, nullable=False),
        sa.Column('LHOOD', sapg.DOUBLE_PRECISION, nullable=False),
        sa.Column('KAPPA_MLE', sapg.DOUBLE_PRECISION, nullable=False),
        sa.Column('OMEGA_FWR_MLE', sapg.DOUBLE_PRECISION, nullable=False),
        sa.Column('OMEGA_CDR_MLE', sapg.DOUBLE_PRECISION, nullable=False),
        sa.Column('WRC_2_MLE', sapg.DOUBLE_PRECISION, nullable=False),
        sa.Column('GYW_0_MLE', sapg.DOUBLE_PRECISION, nullable=False),
        sa.Column('WA_1_MLE', sapg.DOUBLE_PRECISION, nullable=False),
        sa.Column('TW_0_MLE', sapg.DOUBLE_PRECISION, nullable=False),
        sa.Column('SYC_2_MLE', sapg.DOUBLE_PRECISION, nullable=False),
        sa.Column('GRS_0_MLE', sapg.DOUBLE_PRECISION, nullable=False),
        sa.Column('TREE', sa.String(), nullable=False),
    )

    return SimpleNamespace(
        antigen=antigen,
        tree2=tree2,
    )


def load_data_tree2_AB_results(conn: sae.Connection, antigen_df: pd.DataFrame, tree_dir: str, dst: sa.Table):
    # fdf: pd.DataFrame = pd.read_csv(tree_f, sep='\t', skiprows=[1, ],
    #                                 dtype={'CLONE': np.int32, })
    # sdf: pd.DataFrame = pd.read_sql_table(table_name=dst.name, schema=dst.schema, con=conn)
    # fdf[~(fdf.CLONE.isin(sdf.CLONE)) & (fdf['CLONE'] != 'REPERTOIRE')] \
    #     .to_sql(name=dst.name, schema=dst.schema, con=conn, index=False, if_exists='append')
    for (_, ag_row) in antigen_df.iterrows():
        antigen_id: int = ag_row['id']
        antigen: str = ag_row['antigen']

        antigen_tree_dir = os.path.join(tree_dir, 'AB_AG_data_{0}'.format(antigen));
        if os.path.isdir(antigen_tree_dir):
            load_data_tree2_antigen(conn, antigen_id, antigen, antigen_tree_dir, dst)
        else:
            sys.stderr.write(
                "Antigen '{ag}' tree data dir '{dir}' not found.\n".format(ag=antigen, dir=antigen_tree_dir))
    k = 11


def load_data_tree2_antigen(conn: sae.Connection, antigen_id: int, antigen: str, antigen_tree_dir: str, dst: sa.Table):
    try:
        tree_csv_fn: str = os.path.join(
            antigen_tree_dir,
            'AB_AG_data_{0}_db-pass_parse-select_clone-pass_germ-pass_igphyml-pass.tab'.format(antigen))
        fdf: pd.DataFrame = pd.read_csv(tree_csv_fn, delimiter='\t')
        fdf['antigen_id'] = [antigen_id] * len(fdf)
        fdf['antigen'] = [antigen] * len(fdf)
        fdf.to_sql(name=dst.name, schema=dst.schema, con=conn, index=False, if_exists='append')
        fdf_tree2: pd.DataFrame = fdf['CLONE']
    except Exception:
        # ex_type, ex_value, ex_trace = sys.exc_info()
        sys.stderr.write(traceback.format_exc())
        raise  # ex_type(ex_value).with_traceback(ex_trace)


@click.group(cls=DefaultGroup, default='main')
def cli():
    pass


@cli.command()
@click.pass_context
@click.option('--conn-str', 'conn_str',
              help="Database connection string .",
              type=click.STRING)
@click.option('--tree-dir', 'tree_dir',
              help='Directory with igphyml results',
              type=click.Path(exists=True, readable=True, dir_okay=True))
def main(ctx, conn_str, tree_dir: str):
    s = build_schema()
    engine = sa.create_engine(conn_str)
    meta_mlb.create_all(bind=engine)

    with engine.begin() as conn:
        antigen_sdf: pd.DataFrame = pd.read_sql_table(table_name=s.antigen.name, schema=s.antigen.schema, con=conn)

        load_data_tree2_AB_results(conn, antigen_sdf, tree_dir, s.tree2)

    k = 11


if __name__ == '__main__':
    cli()
