#!/usr/bin/env python

import os.path
from datetime import datetime
from types import SimpleNamespace

import numpy as np
import pandas as pd
import re

import click
from click_default_group import DefaultGroup

import sqlalchemy as sa
import sqlalchemy.dialects.postgresql as sapg
import sqlalchemy.engine as sae
import sqlalchemy.orm as sao

meta_public = sa.MetaData(schema='public')
meta_v2 = sa.MetaData(schema='db_v2')
Base_v2 = sao.declarative_base(metadata=meta_v2)

scheme_list = ['chothia', 'imgt', 'kabat']
chain_list = ['heavy', 'light']


class NumberingSchemeMapper:
    def __init__(self, numbering_scheme_f):
        # We have to replace column names, due to some confusion in the file headers
        self._ns_df: pd.DataFrame = pd.read_csv(
            numbering_scheme_f, sep='\t', skiprows=[0, 1, ], dtype=str,
            names=['aho', 'chothia', 'chothia_enh', 'kabat', 'imgt'])

    def map_position(self, src_scheme: str, src_pos_name: str, pos_shift: int, tgt_scheme: str):
        src_pos_list: list[str] = self._ns_df.loc[:, src_scheme].to_list()
        pos_idx = src_pos_list.index(src_pos_name)
        # if pos_idx == -1:
        #     raise ValueError(f"Position '{src_pos_name}' not found in numbering scheme definitions.")
        tgt_scheme_col_idx = self._ns_df.columns.to_list().index(tgt_scheme)
        # if tgt_scheme_col_idx == -1:
        #     raise ValueError(f"Target scheme '{tgt_scheme}' not fount in numbering scheme definitions.")
        return self._ns_df.iloc[pos_idx + pos_shift, tgt_scheme_col_idx]

    def map_first(self, scheme_name: str, chain_sel: str):
        pos_name_res: str = self._ns_df.loc[self._ns_df[scheme_name].str.startswith(chain_sel), scheme_name].iloc[0]
        return pos_name_res

    def map_last(self, scheme_name: str, chain_sel: str):
        pos_name_res: str = self._ns_df.loc[self._ns_df[scheme_name].str.startswith(chain_sel), scheme_name].iloc[-1]
        return pos_name_res

    @classmethod
    def simplify_position_name(cls, pos_name):
        """
        :param pos_name: position name from (.tsv) text files with
                         numbering scheme definitions/mapping and cdr definitions
        :return: simplified position name
        """
        pos_m: re.Match = re.search(r'[LH]:(\d+):([~A-Z])', pos_name)
        pos_res: str = pos_m.group(1) if pos_m.group(2) == '~' \
            else f'{pos_m.group(1)}{pos_m.group(2)}'
        return pos_res


class VdRegion:
    def __init__(self, chain: str, type: str, name: str, order: int,
                 pos_start: str, pos_start_shift: int, pos_end: str, pos_end_shift: int):
        self.chain = chain
        self.type = type
        self.name = name
        self.order = order
        self.pos_start = pos_start
        self.pos_start_shift = pos_start_shift
        self.pos_end = pos_end
        self.pos_end_shift = pos_end_shift

    def __repr__(self):
        return f"chain: '{self.chain}', type: '{self.type}', name: '{self.name}', order: {self.order}, " \
               f"start: '{self.pos_start}' {self.pos_start_shift}, end: '{self.pos_end}' {self.pos_end_shift}"

    def __str__(self):
        return self.__repr__()


class CdrLayoutBuilder:
    """
    Generalized description of the area of regions in terms
    of positions names as defined in file cdr_definitions.tsv
    """
    cdr_region_list = [
        VdRegion('Light', 'framework', 'FR1', 1, 'start', +1, 'L1 cdr_start', -1),
        VdRegion('Heavy', 'framework', 'FR1', 1, 'start', +1, 'H1 cdr_start', -1),

        VdRegion('Light', 'cdr', 'CDR1', 2, 'L1 cdr_start', 0, 'L1 cdr_end', 0),
        VdRegion('Heavy', 'cdr', 'CDR1', 2, 'H1 cdr_start', 0, 'H1 cdr_end', 0),

        VdRegion('Light', 'framework', 'FR2', 3, 'L1 cdr_end', +1, 'L2 cdr_start', -1),
        VdRegion('Heavy', 'framework', 'FR2', 3, 'H1 cdr_end', +1, 'H2 cdr_start', -1),

        VdRegion('Light', 'cdr', 'CDR2', 4, 'L2 cdr_start', 0, 'L2 cdr_end', 0),
        VdRegion('Heavy', 'cdr', 'CDR2', 4, 'H2 cdr_start', 0, 'H2 cdr_end', 0),

        VdRegion('Light', 'framework', 'FR3', 5, 'L2 cdr_end', +1, 'L3 cdr_start', -1),
        VdRegion('Heavy', 'framework', 'FR3', 5, 'H2 cdr_end', +1, 'H3 cdr_start', -1),

        VdRegion('Light', 'cdr', 'CDR3', 6, 'L3 cdr_start', 0, 'L3 cdr_end', 0),
        VdRegion('Heavy', 'cdr', 'CDR3', 6, 'H3 cdr_start', 0, 'H3 cdr_end', 0),

        VdRegion('Light', 'framework', 'FR4', 7, 'L3 cdr_end', +1, 'end', -1),
        VdRegion('Heavy', 'framework', 'FR4', 7, 'H3 cdr_end', +1, 'end', -1),
    ]

    def __init__(self, cdr_f, numbering_scheme_f):
        # We have to replace column names, due to some confusion in the file headers
        self._cdr_df: pd.DataFrame = pd.read_csv(
            cdr_f, sep='\t', skiprows=[0, 1, 2, ], dtype=str,
            names=['chothia:chothia', 'aroop:chothia', 'north:aho', 'martin:aho', 'kabat:aho'])
        self._scheme_mapper = NumberingSchemeMapper(numbering_scheme_f)

    def get_cdr_region_list(self, scheme_name: str, cdr_name: str, scheme_layout_table: sa.Table):
        """
        List of regions for numbering scheme and cdr definition arguments.
        :param scheme_name:
        :param cdr_name:
        :param scheme_layout_table: to generate Insert objects
        :return: List of regions (preferable records of table 'scheme_layout')
        """

        def get_position(region_pos_name: str, region_pos_shift: int) -> str:
            # search region_pos_name in self._cdr_df row names
            cdr_row_idx = self._cdr_df.index.to_list().index(region_pos_name)
            # if cdr_row_idx == -1:
            #     raise ValueError(f"Region position name '{region_pos_name}' not found in cdr definitions.")

            tgt_pos_name: str
            if f'{cdr_name}:{scheme_name}' in self._cdr_df.columns:
                tgt_pos_name = self._cdr_df.loc[region_pos_name, f'{cdr_name}:{scheme_name}']
            else:
                cdr_col_idx = [cdr_col_name.split(':')[0] for cdr_col_name in self._cdr_df.columns].index(cdr_name)
                # if cdr_col_idx == -1:
                #     raise ValueError(f"CDR '{cdr_name}' not found in cdr definitions.")

                cdr_src_scheme = self._cdr_df.columns[cdr_col_idx].split(':')[1]
                cdr_pos_name = self._cdr_df.iloc[cdr_row_idx, cdr_col_idx]
                tgt_pos_name = self._scheme_mapper \
                    .map_position(cdr_src_scheme, cdr_pos_name, region_pos_shift, scheme_name)

            return tgt_pos_name

        def build_record(region: VdRegion):
            pos_start = NumberingSchemeMapper.simplify_position_name(
                self._scheme_mapper.map_first(scheme_name, region.chain[0]) if region.pos_start == 'start'
                else get_position(region.pos_start, region.pos_start_shift))

            pos_end = NumberingSchemeMapper.simplify_position_name(
                self._scheme_mapper.map_last(scheme_name, region.chain[0]) if region.pos_end == 'end'
                else get_position(region.pos_end, region.pos_end_shift))

            rec_res = scheme_layout_table.insert().values(
                scheme=scheme_name,
                cdr=cdr_name,
                type=region.type,
                name=region.name,
                chain=region.chain,
                order=region.order,
                position_start_name=pos_start,
                position_end_name=pos_end
            )
            return rec_res

        rec_list_res = [build_record(region) for region in CdrLayoutBuilder.cdr_region_list]
        return rec_list_res


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
    list_version = sa.Table(
        'list_version', meta_v2,
        sa.Column('list_id', sa.Integer, primary_key=True, autoincrement=True),
        sa.Column('name', sa.String(), nullable=False, unique=True),
        sa.Column('version', sa.String(), nullable=False),
        # TODO: Use timestamp column and update it with triggers on list's tables
    )

    antigen = sa.Table(
        'antigen', meta_v2,
        sa.Column('antigen_id', sa.Integer, primary_key=True, autoincrement=True),
        sa.Column('antigen', sa.String(), nullable=False, unique=True),
        sa.Column('antigen_ncbi_id', sa.String(), nullable=True),
        sa.Column('antigen_gene_symbol', sa.String(), nullable=True)
    )

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

    antibody_anarci_tables: dict[str, dict[str, any]] = \
        dict([(
            scheme,
            dict([(
                chain,
                build_antibody_anarci_table(antibody, scheme, chain, anarci_path)
            )  # ANARCI results file name
                for chain in chain_list])
        ) for scheme in scheme_list])

    scheme = sa.Table(
        'scheme', meta_v2,
        sa.Column('scheme_id', sa.Integer, primary_key=True, autoincrement=True),
        sa.Column('scheme', sa.String(), nullable=False, unique=True)
    )

    cdr = sa.Table(
        'cdr', meta_v2,
        sa.Column('cdr_id', sa.Integer, primary_key=True, autoincrement=True),
        sa.Column('cdr', sa.String(), nullable=False, unique=True),
    )

    scheme_layout = sa.Table(
        'scheme_layout', meta_v2,
        sa.Column('region_id', sa.Integer, primary_key=True, autoincrement=True),
        sa.Column('scheme', sa.String, sa.ForeignKey(scheme.columns['scheme'])),
        sa.Column('cdr', sa.String, sa.ForeignKey(cdr.columns['cdr'])),
        sa.Column('type', sa.String(), nullable=False),
        sa.Column('name', sa.String(), nullable=False),
        sa.Column('chain', sa.String(), nullable=False),
        sa.Column('order', sa.Integer, nullable=False),
        sa.Column('position_start_name', sa.String(), nullable=True),
        sa.Column('position_end_name', sa.String(), nullable=True),
    )

    return SimpleNamespace(
        list_version=list_version,
        antigen=antigen,
        antibody=antibody,
        antibody2antigen=antibody2antigen,
        tree=tree,
        antibody2tree=antibody2tree,
        antibody_anarci_tables=antibody_anarci_tables,
        scheme=scheme,
        cdr=cdr,
        scheme_layout=scheme_layout,
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
    ids = set(vr_re.findall(tree_txt))
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


def load_data_scheme(conn: sae.Connection, scheme: sa.Table):
    ins_dict = {
        'imgt': scheme.insert().values(scheme='imgt'),
        'aho': scheme.insert().values(scheme='aho'),
        'chothia': scheme.insert().values(scheme='chothia'),
        'kabat': scheme.insert().values(scheme='kabat'),
    }
    for (scheme_name, scheme_ins) in ins_dict.items():
        scheme_check_stmt = sa.select(scheme).filter_by(scheme=scheme_name)
        scheme_check_res = conn.execute(scheme_check_stmt).one_or_none()
        if not scheme_check_res:
            conn.execute(scheme_ins)


def load_data_cdr(conn: sae.Connection, cdr: sa.Table):
    ins_dict = {
        'chothia': cdr.insert().values(cdr='chothia'),
        'aroop': cdr.insert().values(cdr='aroop'),
        'north': cdr.insert().values(cdr='north'),
        'martin': cdr.insert().values(cdr='martin'),
        'kabat': cdr.insert().values(cdr='kabat'),
    }
    for (cdr_name, cdr_ins) in ins_dict.items():
        cdr_check_stmt = sa.select(cdr).filter_by(cdr=cdr_name)
        cdr_check_res = conn.execute(cdr_check_stmt).one_or_none()
        if not cdr_check_res:
            conn.execute(cdr_ins)


def load_data_scheme_layout_for_scheme_cdr(
        conn: sae.Connection, cdr_layout_builder: CdrLayoutBuilder,
        scheme_layout_table: sa.Table, scheme_name: str, cdr_name: str) -> None:
    rec_list = cdr_layout_builder.get_cdr_region_list(scheme_name, cdr_name, scheme_layout_table)
    for rec in rec_list:
        conn.execute(rec)


def load_data_scheme_layout(conn: sae.Connection,
                            scheme_table: sa.Table, cdr_table: sa.Table,
                            cdr_layout_builder: CdrLayoutBuilder, scheme_layout_table: sa.Table):
    for (scheme_id, scheme_name) in conn.execute(sa.select(scheme_table)).all():
        for (crd_id, cdr_name) in conn.execute(sa.select(cdr_table)).all():
            layout_check_res = conn.execute(
                sa.select(scheme_layout_table).filter_by(scheme=scheme_name, cdr=cdr_name)).all()
            if len(layout_check_res) == 0:
                load_data_scheme_layout_for_scheme_cdr(
                    conn, cdr_layout_builder, scheme_layout_table, scheme_name, cdr_name)
            else:
                print(f"There are already {len(layout_check_res)} layout records "
                      f"for scheme '{scheme_name}' cdr '{cdr_name}'.")


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
@click.pass_context
@click.option('--conn-str', 'conn_str',
              help="Database connection string .",
              type=click.STRING)
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
@click.option('--numbering-scheme', 'numbering_scheme_f',
              help='Text (.tsv) file with mapping numbering schemes positions',
              type=click.File('r'))
@click.option('--cdr', 'cdr_f',
              help='Text (.tsv) file with cdr regions definitions',
              type=click.File('r'))
def main(ctx, conn_str, ag_f, ab_f, ab2ag_f, tree_f, anarci_path, numbering_scheme_f, cdr_f):
    s = build_schema(anarci_path)

    engine = sa.create_engine(conn_str)

    meta_v2.create_all(bind=engine)

    with engine.begin() as conn:
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

        load_data_antibody_anarci(conn, antibody_fdf, s.antibody_anarci_tables)

        load_data_tree(conn, antibody_fdf, tree_f, s.tree)
        tree_sdf: pd.DataFrame = pd.read_sql_table(table_name=s.tree.name, schema=s.tree.schema, con=conn)
        load_data_ab2tr(conn, antibody_fdf, tree_sdf, s.antibody2tree)

        load_data_scheme(conn, s.scheme)
        load_data_cdr(conn, s.cdr)
        cdr_layout_builder = CdrLayoutBuilder(cdr_f, numbering_scheme_f)
        load_data_scheme_layout(conn, s.scheme, s.cdr, cdr_layout_builder, s.scheme_layout, )

        k = 11


if __name__ == '__main__':
    cli()
