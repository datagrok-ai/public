#!/usr/bin/env python
import itertools
import os.path
import io
import sys
from datetime import datetime
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd
import re

import click
from click_default_group import DefaultGroup

import sqlalchemy as sa
import sqlalchemy.dialects.postgresql as sapg
import sqlalchemy.engine as sae
import sqlalchemy.orm as sao

import Bio.Phylo as phy
import Bio.Phylo.NewickIO as nwkio
import Bio.Phylo.Newick as nwk

meta_public = sa.MetaData(schema='public')
meta_db_v2 = sa.MetaData(schema='db_v2')
meta_mlb = sa.MetaData(schema='mlb')


def build_schema():
    antigen = sa.Table(
        'antigen', meta_mlb,
        sa.Column('antigen_id', sa.Integer, primary_key=True, autoincrement=True),
        sa.Column('antigen', sa.String(), nullable=False, unique=True),
        sa.Column('antigen_ncbi_id', sa.String(), nullable=True),
        sa.Column('antigen_gene_symbol', sa.String(), nullable=True)
    )

    antibody = sa.Table(  # MLB
        'mlb_main', meta_mlb,
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
        'antibody2antigen', meta_mlb,
        sa.Column('ab2ag_id', sa.Integer, primary_key=True, autoincrement=True),
        sa.Column('antigen', sa.String(), sa.ForeignKey(antigen.columns['antigen'])),
        sa.Column('v_id', sa.String(), sa.ForeignKey(antibody.columns['v_id'])),
    )

    tree_col_list: list[sa.Column] = [
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
    ]

    tree = sa.Table(
        'tree', meta_mlb,
        *[col._copy() for col in tree_col_list]
    )

    # obsolete
    # helper table to search tree/clone by antigen
    # antibody2tree = sa.Table(
    #     'antibody2tree', meta_mlb,
    #     sa.Column('ab2tr_id', sa.Integer, primary_key=True, autoincrement=True),
    #     sa.Column('CLONE', sa.Integer, sa.ForeignKey(tree.columns['CLONE'])),
    #     sa.Column('v_id', sa.String(), sa.ForeignKey(antibody.columns['v_id'])),
    # )
    tree2 = sa.Table(
        'tree2', meta_mlb,
        sa.Column('tree_id', sa.Integer, primary_key=True),
        sa.Column('antigen_id', sa.Integer, nullable=False),
        sa.Column('antigen', sa.String(), nullable=False),
        *[col._copy() for col in tree_col_list]
    )

    return SimpleNamespace(
        antigen=antigen,
        antibody=antibody,
        antibody2antigen=antibody2antigen,
        tree2=tree2
    )


def check_antigen(
        ag_name: str, ab2ag: pd.DataFrame, tree2: pd.DataFrame, all_vr_set: set[str], res_col_list: list[pd.Series]
) -> (int, pd.DataFrame):
    ag_vr_set: set[str] = set()  # all VRs associated with the antigen (by ab2ag table)
    ag_ab2ag_df: pd.DataFrame = ab2ag[ab2ag.antigen == ag_name]
    for (_, ab2ag_row) in ag_ab2ag_df.iterrows():
        v_id = ab2ag_row['v_id']
        ag_vr_set.add(v_id)

    tree_vr_dict: {} = {}  # all VRs associated with the antigen (by tree2 table TREE column parsed)

    def get_leaf_list(clade: nwk.Clade) -> list[str]:
        if clade.clades and len(clade.clades) > 0:
            # TODO: Merge generator of generators
            # yield from itertools.chain.from_iterable([get_leaf_list(child_clade) for child_clade in clade.clades])
            # res = [j for child_clade in clade.clades for j in i]
            # # [j for i in x for j in i]
            # return res
            for child_clade in clade.clades:
                yield from get_leaf_list(child_clade)
        else:  # leaf
            yield clade

    def get_vid(leaf_name: str) -> str:
        leaf_name_part_list = leaf_name.split('|')
        if len(leaf_name_part_list) != 4:
            raise ValueError("Tree / clone leaf name '{0}' is of invalid format.".format(leaf_name))
        return leaf_name_part_list[2]

    ag_vid_set: set = set()
    ag_res_list: list[dict[str, Any]] = []
    ag_tree2_df: pd.DataFrame = tree2[tree2.antigen == ag_name]
    for (_, tree2_row) in ag_tree2_df.iterrows():
        tree_clone = tree2_row.CLONE
        if tree_clone != 'REPERTOIRE':
            tree_nwk = tree2_row.TREE
            tree_obj = phy.read(io.StringIO(tree_nwk), 'newick')
            tree_leaf_gen = get_leaf_list(tree_obj.clade)
            tree_vid_gen = (get_vid(tree_leaf.name) for tree_leaf in tree_leaf_gen
                            if not tree_leaf.name.endswith('_GERM'))
            tree_vid_list = list(tree_vid_gen)
            for v_id in tree_vid_list:
                if v_id not in ag_vid_set:  # prevent duplicate VR in res
                    ag_vid_set.add(v_id)
                    if v_id not in ag_vr_set:
                        ag_res_list.append({
                            'antigen': ag_name,
                            'v_id': v_id,  # v_id from parsed clone tree
                            'in_antigen': '1' if v_id in ag_vr_set else '0',
                            'in_all': '1' if v_id in all_vr_set else '0',
                            'CLONE': tree_clone,
                            'TREE': tree_nwk,
                        })
    # ag_res_df = pd.DataFrame(ag_res_list, columns=res_col_list)
    ag_res_df = pd.DataFrame(data=dict([
        (col.name, pd.Series([row[col.name] for row in ag_res_list], dtype=col.dtype)) for col in res_col_list
    ]))
    return (len(ag_vid_set), ag_res_df)


def read_sql_table(table: sa.Table, con: Any) -> pd.DataFrame:
    return pd.read_sql_table(table_name=table.name, schema=table.schema, con=con)


@click.group(cls=DefaultGroup, default='main')
def cli():
    pass


@cli.command()
@click.pass_context
@click.option('--conn-str', 'conn_str',
              help="Database connection string",
              type=click.STRING)
@click.option('-o', '--out', 'out_f',
              help='Output report CSV file name',
              type=click.File('w'))
def main(ctx, conn_str, out_f):
    s = build_schema()
    engine = sa.create_engine(conn_str)

    # Do not modify the database
    # meta_mlb.create_all(bind=engine)

    with engine.begin() as conn:
        ag_sdf: pd.DataFrame = read_sql_table(s.antigen, con=conn)
        ab_sdf: pd.DataFrame = read_sql_table(s.antibody, con=conn)
        ab2ag_sdf: pd.DataFrame = read_sql_table(s.antibody2antigen, con=conn)
        tree2_sdf: pd.DataFrame = read_sql_table(s.tree2, con=conn)

        res_col_list: list[pd.Series] = [
            pd.Series(name='antigen', dtype=pd.StringDtype.name),
            pd.Series(name='v_id', dtype=pd.StringDtype.name),  # v_id from parsed clone tree
            pd.Series(name='in_antigen', dtype=pd.UInt8Dtype.name),
            pd.Series(name='in_all', dtype=pd.UInt8Dtype.name),
            pd.Series(name='CLONE', dtype=pd.StringDtype.name),
            pd.Series(name='TREE', dtype=pd.StringDtype.name),
        ]
        res_df: pd.DataFrame = pd.DataFrame(
            data=dict([(col.name, pd.Series([], dtype=col.dtype)) for col in res_col_list]))
        all_vr_set: set[str] = set(ab_sdf.loc[:, 'v_id'])
        for (_, ag_row) in ag_sdf.iterrows():
            ag_name: str = ag_row['antigen']

            (ag_vid_count, ag_res_df) = \
                check_antigen(ag_name, ab2ag_sdf, tree2_sdf, all_vr_set, res_col_list)
            sys.stderr.write('')

            if len(ag_res_df) > 0:
                missed_ag_count = len(ag_res_df[ag_res_df['in_antigen'] == 0])
                missed_all_count = len(ag_res_df[ag_res_df['in_all'] == 0])
                sys.stderr.write(
                    ("Antigen '{0}': total {total}, " +
                     "missed with antigen {missed_ag_count}, missed at all {missed_all_count}\n")
                    .format(ag_name, total=ag_vid_count,
                            missed_ag_count=missed_ag_count, missed_all_count=missed_all_count))
            res_df = pd.concat([res_df, ag_res_df], axis=0, ignore_index=True)
            k = 11

        # res_fn = 'check_vr_for_tree.csv'
        res_df.to_csv(out_f, index=False, sep=',', line_terminator='\n')


if __name__ == '__main__':
    cli()
