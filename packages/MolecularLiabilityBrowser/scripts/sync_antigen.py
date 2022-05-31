#!/usr/bin/env python
import os.path
import sys
from types import SimpleNamespace

import pandas as pd
import numpy as np
import re
import orjson as json

import psycopg2 as pg

import click
from click_default_group import DefaultGroup

import sqlalchemy as sa
import sqlalchemy.engine as sae
import sqlalchemy.orm as sao
import sqlalchemy.dialects.postgresql as sapg

meta_public = sa.MetaData(schema='public')
meta_v2 = sa.MetaData(schema='db_v2')
Base_v2 = sao.declarative_base(metadata=meta_v2)

scheme_list = ['chothia', 'imgt', 'kabat']
chain_list = ['heavy', 'light']


def build_antibody_anarci_table(antibody: sa.Table, scheme: str, chain: str, anarci_path: str) -> sa.Table:
    fn = os.path.join(anarci_path, 'anarci_{0}_{1}_{2}.csv'.format(scheme, chain, {'heavy': 'H', 'light': 'KL'}[chain]))

    table_name = 'antibody_anarci_{0}_{1}'.format(scheme, chain)
    # Dynamic schema in SQLAlchemy
    # https://sparrigan.github.io/sql/sqla/2016/01/03/dynamic-tables.html

    table_columns = [
        sa.Column('Id', sa.String(), sa.ForeignKey(antibody.columns['v_id']), primary_key=True),
        sa.Column('domain_no', sa.Integer, nullable=False),
        sa.Column('hmm_species', sa.String(), nullable=False),
        sa.Column('chain_type', sa.String(), nullable=False),
        sa.Column('e-value', sapg.DOUBLE_PRECISION, nullable=False),
        sa.Column('score', sapg.DOUBLE_PRECISION, nullable=False),
        sa.Column('seqstart_index', sa.Integer, nullable=False),
        sa.Column('seqend_index', sa.Integer, nullable=False),
        sa.Column('identity_species', sa.String(), nullable=True),
        sa.Column('v_gene', sa.String(), nullable=True),
        sa.Column('v_identity', sa.Float, nullable=False),
        sa.Column('j_gene', sa.String(), nullable=True),
        sa.Column('j_identity', sa.Float, nullable=False),
    ]

    fdf: pd.DataFrame = pd.read_csv(fn, dtype=str)
    fdf_cols = list(fdf)  # list columns' names

    first_seq_col = fdf_cols.index('1')  # Raises ValueError if the value '1' is not present
    # add position columns to sqlalchemy table
    for col_i in range(first_seq_col, len(fdf_cols)):
        col_name: str = fdf_cols[col_i]
        table_columns.append(sa.Column(col_name, sa.String(), default='-', nullable=False))
    table = sa.Table(table_name, meta_v2, *table_columns)

    return SimpleNamespace(src=fdf, dst=table)


def build_schema(anarci_path: str):
    antigen = sa.Table(
        'antigen', meta_v2,
        sa.Column('antigen_id', sa.Integer, primary_key=True, autoincrement=True),
        sa.Column('antigen', sa.String(), nullable=False, unique=True),
        sa.Column('antigen_ncbi_id', sa.String(), nullable=True),
        sa.Column('antigen_gene_symbol', sa.String(), nullable=True))

    antibody = sa.Table(  # MLB
        'mlb_main', meta_public,
        sa.Column('v_id', sa.String(), primary_key=True, nullable=False),
        sa.Column('gdb_id_mappings', sa.String(), nullable=True),
        sa.Column('cdr_length', sapg.DOUBLE_PRECISION, nullable=True),
        sa.Column('surface_cdr_hydrophobicity', sapg.DOUBLE_PRECISION, nullable=True),
        sa.Column('positive_cdr_charge', sapg.DOUBLE_PRECISION, nullable=True),
        sa.Column('negative_cdr_charge', sapg.DOUBLE_PRECISION, nullable=True),
        sa.Column('SFvCSP', sapg.DOUBLE_PRECISION, nullable=True),
        # sa.Column('H Hydroxylysine', sapg.REAL, nullable=True),
        # sa.Column('H Methyllysine', sapg.REAL, nullable=True),
        # sa.Column('H N6-acetyllysine', sapg.REAL, nullable=True),
        # sa.Column('H O-linked_glycosylation', sapg.REAL, nullable=True),
        # sa.Column('H Phosphoserine_Phosphothreonine', sapg.REAL, nullable=True),
        # sa.Column('H Phosphotyrosine', sapg.REAL, nullable=True),
        # sa.Column('H Pyrrolidone_carboxylic_acid', sapg.REAL, nullable=True),
        # sa.Column('H S-palmitoyl_cysteine', sapg.REAL, nullable=True),
        # sa.Column('H SUMOylation', sapg.REAL, nullable=True),
        # sa.Column('H Ubiquitination', sapg.REAL, nullable=True),
        # sa.Column('L Hydroxylysine', sapg.REAL, nullable=True),
        # sa.Column('L Hydroxyproline', sapg.REAL, nullable=True),
        # sa.Column('L Methyllysine', sapg.REAL, nullable=True),
        # sa.Column('L N6-acetyllysine', sapg.REAL, nullable=True),
        # sa.Column('L O-linked_glycosylation', sapg.REAL, nullable=True),
        # sa.Column('L Phosphoserine_Phosphothreonine', sapg.REAL, nullable=True),
        # sa.Column('L Phosphotyrosine', sapg.REAL, nullable=True),
        # sa.Column('L Pyrrolidone_carboxylic_acid', sapg.REAL, nullable=True),
        # sa.Column('L SUMOylation', sapg.REAL, nullable=True),
        # sa.Column('L Ubiquitination', sapg.REAL, nullable=True),
    )

    antibody2antigen = sa.Table(
        'antibody2antigen', meta_v2,
        sa.Column('ab2ag_id', sa.Integer, primary_key=True, autoincrement=True),
        sa.Column('antigen', sa.String(), sa.ForeignKey(antigen.columns['antigen'])),
        sa.Column('v_id', sa.String(), sa.ForeignKey(antibody.columns['v_id'])),
    )

    tree = sa.Table(
        'tree', meta_v2,
        sa.Column('CLONE', sa.Integer, primary_key=True),
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

    # helper table to search tree/clone by antigen
    antibody2tree = sa.Table(
        'antibody2tree', meta_v2,
        sa.Column('ab2tr_id', sa.Integer, primary_key=True, autoincrement=True),
        sa.Column('CLONE', sa.Integer, sa.ForeignKey(tree.columns['CLONE'])),
        sa.Column('v_id', sa.String(), sa.ForeignKey(antibody.columns['v_id'])),
    )

    # antibody_anarci_tables: dict[str, dict[str, any]] = \
    #     dict([(
    #         scheme,
    #         dict([(
    #             chain,
    #             build_antibody_anarci_table(antibody, scheme, chain, anarci_path)
    #         )  # ANARCI results file name
    #             for chain in chain_list])
    #     ) for scheme in scheme_list])

    return SimpleNamespace(
        antigen=antigen,
        antibody=antibody,
        antibody2antigen=antibody2antigen,
        tree=tree,
        antibody2tree=antibody2tree,
        # antibody_anarci_tables=antibody_anarci_tables
    )


def load_data_antibody_anarci(conn: sae.Connection, antibody_df: pd.DataFrame, tables: dict[str, dict[str, any]]):
    """
    There is no place to put position names in database table (only as a part of column name - bad idea).
    So we will not precalculate joined 'seq' column values.
    :param conn: Database connection
    :param antibody_df: antibodies allowed in database
    :param tables: Dictionary with ANARCI data and tables' structures
    """
    for scheme in tables.keys():
        for chain in tables[scheme].keys():
            dst = tables[scheme][chain].dst
            fdf: pd.DataFrame = tables[scheme][chain].src
            sdf: pd.Dataframe = pd.read_sql_table(table_name=dst.name, schema=dst.schema, con=conn)
            fdf[~(fdf.Id.isin(sdf.Id)) & fdf.Id.isin(antibody_df.v_id)] \
                .to_sql(name=dst.name, schema=dst.schema, con=conn, index=False, if_exists='append')


def load_data_tree(conn: sae.Connection, antibody_df: pd.DataFrame, tree_f, dst: sa.Table):
    fdf: pd.DataFrame = pd.read_csv(tree_f, sep='\t', skiprows=[1, ],
                                    dtype={'CLONE': np.int32, })
    sdf: pd.DataFrame = pd.read_sql_table(table_name=dst.name, schema=dst.schema, con=conn)
    fdf[~(fdf.CLONE.isin(sdf.CLONE)) & (fdf['CLONE'] != 'REPERTOIRE')] \
        .to_sql(name=dst.name, schema=dst.schema, con=conn, index=False, if_exists='append')


vr_re = re.compile(r'VR\d+')


def get_v_id_list(tree_txt) -> list[str]:
    ids = vr_re.findall(tree_txt)
    return ids


def load_data_ab2tr(conn: sae.Connection, antibody_df: pd.DataFrame, tree_df, dst: sa.Table):
    ab2tr_records: [{}] = []
    for (_, tree_row) in tree_df.iterrows():
        v_id_list: list[str] = get_v_id_list(tree_row['TREE'])
        for v_id in v_id_list:
            ab2tr_records.append({'CLONE': tree_row['CLONE'], 'v_id': v_id})
    ab2tr_src = pd.DataFrame.from_records(ab2tr_records)

    ab2tr_sdf = pd.read_sql_table(table_name=dst.name, schema=dst.schema, con=conn)
    ab2tr_src[
        ~(ab2tr_src.v_id.isin(ab2tr_sdf.v_id) & ab2tr_src.CLONE.isin(ab2tr_sdf.CLONE)) &
        ab2tr_src.v_id.isin(antibody_df.v_id)
        ] \
        .to_sql(name=dst.name, schema=dst.schema, con=conn, index=False, if_exists='append')
    k = 11


@click.group(cls=DefaultGroup, default='main')
def cli():
    pass


# def ensure_table_exists_antigen(engine: sa.Engine):
#     Antigen.__table__.create(bind=engine, check_first=True);
#
#     # check_res = conn.execute(text(sql_antigen_table_check)).all()
#     # if len(check_res):
#     #     create_res = conn.execute(text(sql_antigen_table_create))
#
#     k = 11


@cli.command()
@click.option('--conn-str', 'conn_str',
              help="Database connection string .",
              type=str)
@click.option('--ag', 'ag_f',
              help='CSV file with antigen list',
              type=click.File('r'))
@click.option('--ab', 'ab_f',
              help='CSV file with antibodies (mlb_main)')
@click.option('--ab2ag', 'ab2ag_f',
              help='CSV file with antibody to antigen map',
              type=click.File('r'))
@click.option('--tree', 'tree_f',
              help='CSV file with trees of antibodies',
              type=click.File('r'))
@click.option('--anarci-path', 'anarci_path',
              help='Path to MA ANARCI .csv files',
              type=click.Path(exists=True))
@click.pass_context
def main(ctx, conn_str, ag_f, ab_f, ab2ag_f, tree_f, anarci_path):
    s = build_schema(anarci_path)

    engine = sa.create_engine(conn_str)

    meta_v2.create_all(bind=engine)

    with engine.connect() as conn:
        antigen_fdf: pd.DataFrame = pd.read_csv(ag_f, dtype=str)
        antigen_sdf: pd.DataFrame = pd.read_sql_table(table_name=s.antigen.name, schema=s.antigen.schema, con=conn)

        antigen_fdf[~antigen_fdf.antigen.isin(antigen_sdf.antigen)] \
            .to_sql(name=s.antigen.name, schema=s.antigen.schema, con=conn, index=False, if_exists='append',
                    dtype={'antigen': sa.String,
                           'antigen_ncbi_id': sa.String,
                           'antigen_gene_symbol': sa.String})

        antibody_fdf: pd.DataFrame = pd.read_csv(ab_f, dtype=str)
        antibody_fdf.drop(['antigen_list', 'antigen_ncbi_id', 'antigen_gene_symbol', ], axis=1, inplace=True)
        antibody_fdf.rename({'sfvcsp': 'SFvCSP'}, axis=1, inplace=True)
        antibody_sdf: pd.DataFrame = pd.read_sql_table(table_name=s.antibody.name, schema=s.antibody.schema, con=conn)
        antibody_fdf[~antibody_fdf.v_id.isin(antibody_sdf.v_id)] \
            .to_sql(name=s.antibody.name, schema=s.antibody.schema, con=conn, index=False, if_exists='append')

        ab2ag_fdf = pd.read_csv(ab2ag_f)
        ab2ag_fdf.antigen_list = [
            re.sub(r'[{}]', '', v) if v != '[Not Found]' else None
            for v in ab2ag_fdf.antigen_list]

        ab2ag_sdf = pd.read_sql_table(table_name=s.antibody2antigen.name, schema=s.antibody2antigen.schema, con=conn)
        # ab2ag_src = pd.DataFrame(data=None, columns=ab2ag_sdf.columns)

        ab2ag_records: list[{}] = []
        for (_, ab2ag_row) in ab2ag_fdf.iterrows():
            if ab2ag_row['antigen_list']:
                v_id = ab2ag_row['v_id']
                antigen_list_str = ab2ag_row.antigen_list
                antigen_list = [ag.strip() for ag in antigen_list_str.split(',')] \
                    if antigen_list_str is not None else []
                for antigen_name in antigen_list:
                    ab2ag_records.append({'antigen': antigen_name, 'v_id': v_id})

        ab2ag_src = pd.DataFrame.from_records(ab2ag_records)
        ab2ag_src[
            ~(ab2ag_src.v_id.isin(ab2ag_sdf.v_id) & ab2ag_src.antigen.isin(ab2ag_sdf.antigen)) &
            # There are many v_id in files-extra/antibody2antigen.csv that are not in mlb.csv
            ab2ag_src.v_id.isin(antibody_fdf.v_id)
            ] \
            .to_sql(name=s.antibody2antigen.name, schema=s.antibody2antigen.schema, con=conn, index=False,
                    if_exists='append')

        # load_data_antibody_anarci(conn, antibody_fdf, s.antibody_anarci_tables)

        load_data_tree(conn, antibody_fdf, tree_f, s.tree)
        tree_sdf: pd.DataFrame = pd.read_sql_table(table_name=s.tree.name, schema=s.tree.schema, con=conn)
        load_data_ab2tr(conn, antibody_fdf, tree_sdf, s.antibody2tree)

        k = 11


if __name__ == '__main__':
    cli()
