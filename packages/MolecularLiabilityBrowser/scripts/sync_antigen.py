from types import SimpleNamespace

import pandas as pd
import re
import orjson as json

import psycopg2 as pg

import click
from click_default_group import DefaultGroup

import sqlalchemy as sa
import sqlalchemy.orm as sao
import sqlalchemy.dialects.postgresql as sapg

Base = sao.declarative_base()
meta_public = sa.MetaData(schema='public')
meta_v2 = sa.MetaData(schema='db_v2')

Mlb = sa.Table()

antigen = sa.Table(
    'antigen', meta_v2,
    sa.Column('antigen_id', sa.Integer, primary_key=True, autoincrement=True),
    sa.Column('antigen', sa.String(), nullable=False, unique=True),
    sa.Column('antigen_ncbi_id', sa.String(), nullable=True),
    sa.Column('antigen_gene_symbol', sa.String(), nullable=True))

antibody = sa.Table(  # MLB
    'mlb_main', meta_public,
    sa.Column('v_id', sa.String(), primary_key=True, nullable=False),
    sa.Column('ngl', sa.String(), nullable=True),
    sa.Column('gdb_id_mappings', sa.String(), nullable=True),
    sa.Column('cdr_length', sapg.DOUBLE_PRECISION, nullable=True),
    sa.Column('surface_cdr_hydrophobicity', sapg.DOUBLE_PRECISION, nullable=True),
    sa.Column('positive_cdr_charge', sapg.DOUBLE_PRECISION, nullable=True),
    sa.Column('negative_cdr_charge', sapg.DOUBLE_PRECISION, nullable=True),
    sa.Column('sfvcsp', sapg.DOUBLE_PRECISION, nullable=True),
    sa.Column('H Hydroxylysine', sapg.REAL, nullable=True),
    sa.Column('H Methyllysine', sapg.REAL, nullable=True),
    sa.Column('H N6-acetyllysine', sapg.REAL, nullable=True),
    sa.Column('H O-linked_glycosylation', sapg.REAL, nullable=True),
    sa.Column('H Phosphoserine_Phosphothreonine', sapg.REAL, nullable=True),
    sa.Column('H Phosphotyrosine', sapg.REAL, nullable=True),
    sa.Column('H Pyrrolidone_carboxylic_acid', sapg.REAL, nullable=True),
    sa.Column('H S-palmitoyl_cysteine', sapg.REAL, nullable=True),
    sa.Column('H SUMOylation', sapg.REAL, nullable=True),
    sa.Column('H Ubiquitination', sapg.REAL, nullable=True),
    sa.Column('L Hydroxylysine', sapg.REAL, nullable=True),
    sa.Column('L Hydroxyproline', sapg.REAL, nullable=True),
    sa.Column('L Methyllysine', sapg.REAL, nullable=True),
    sa.Column('L N6-acetyllysine', sapg.REAL, nullable=True),
    sa.Column('L O-linked_glycosylation', sapg.REAL, nullable=True),
    sa.Column('L Phosphoserine_Phosphothreonine', sapg.REAL, nullable=True),
    sa.Column('L Phosphotyrosine', sapg.REAL, nullable=True),
    sa.Column('L Pyrrolidone_carboxylic_acid', sapg.REAL, nullable=True),
    sa.Column('L SUMOylation', sapg.REAL, nullable=True),
    sa.Column('L Ubiquitination', sapg.REAL, nullable=True),
)

antibody2antigen = sa.Table(
    'antibody2antigen', meta_v2,
    sa.Column('ab2ag_id', sa.Integer, primary_key=True, autoincrement=True),
    sa.Column('antigen', sa.String(), sa.ForeignKey(antigen.columns['antigen'])),
    sa.Column('v_id', sa.String(), sa.ForeignKey(antibody.columns['v_id'])),
)


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
@click.option('--ag', 'ag_f',
              help='CSV file with antigen list',
              type=click.File('r'))
@click.option('--ab2ag', 'ab2ag_f',
              help='CSV file with antibody to antigen map',
              type=click.File('r'))
@click.pass_context
def main(ctx, ag_f, ab2ag_f):
    # conn: pg.connection = pg.connect({
    #     'host': 'localhost',
    #     'port': '5432',
    #     'user': 'postgres',
    #     'password': 'postgres',
    #     'dbname': 'mlb_db_v2', })

    engine = sa.create_engine('postgresql://postgres:postgres@localhost:5432/mlb_db_v2')
    meta_v2.create_all(bind=engine)

    with engine.connect() as conn:
        antigen_fdf: pd.DataFrame = pd.read_csv(ag_f, dtype=str)
        antigen_sdf: pd.DataFrame = pd.read_sql_table(table_name=antigen.name, schema=antigen.schema, con=conn)

        antigen_fdf[~antigen_fdf.antigen.isin(antigen_sdf.antigen)].to_sql(
            name=antigen.name, schema=antigen.schema, con=conn, index=False, if_exists='append',
            dtype={'antigen': sa.String,
                   'antigen_ncbi_id': sa.String,
                   'antigen_gene_symbol': sa.String})

        ab2ag_fdf = pd.read_csv(ab2ag_f)
        ab2ag_fdf.antigen_list = [
            re.sub(r'[{}]', '', v) if v != '[Not Found]' else None
            for v in ab2ag_fdf.antigen_list]

        ab2ag_sdf = pd.read_sql_table(table_name=antibody2antigen.name, schema=antibody2antigen.schema, con=conn)

        for (_, ab2ag_row) in ab2ag_fdf.iterrows():
            if ab2ag_row['antigen_list']:
                v_id = ab2ag_row['v_id']
                antigen_list_str = ab2ag_row.antigen_list
                antigen_list = [ag.strip() for ag in antigen_list_str.split(',')] \
                    if antigen_list_str is not None else []


        k = 11


if __name__ == '__main__':
    cli()
