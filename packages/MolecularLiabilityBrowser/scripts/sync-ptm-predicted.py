#!/usr/bin/env python

import os.path
from types import SimpleNamespace

import sys
import pandas as pd
import orjson as json
from tqdm import tqdm
from numpy import nan, argmax
from string import ascii_letters

import click
from click_default_group import DefaultGroup

import sqlalchemy as sa
import sqlalchemy.dialects.postgresql as sapg
import sqlalchemy.engine as sae
import sqlalchemy.orm as sao

meta_public = sa.MetaData(schema='public')
meta_mlb = sa.MetaData(schema='mlb')
Base_v2 = sao.declarative_base(metadata=meta_mlb)


def vladimir_script(jsons_sdf: pd.DataFrame, ptm_map: dict[str, str], cdr_map: dict[str, str], dst_sdf: pd.DataFrame):
    def column_ensure(df, colName, value=nan) -> None:
        if colName not in df.columns:
            df[colName] = value

    def key_get_or_create(mapping, key, value):
        if key not in mapping:
            mapping[key] = value
            # sys.stderr.write(f"warn: mapping key '{key}' not found and added")
        return mapping[key]

    def ptm_in_cdr_ranges(ranges, ptm_indexes):
        result = []
        for index, probability in ptm_indexes:
            for low, high in ranges:
                if low <= index <= high:
                    result.append(round(probability, 4) if probability != 100 else 1.0)
        return result

    ptm_mapping = ptm_map.copy()
    cdr_mapping = cdr_map.copy()
    ptm_ascii_index = max((ascii_letters.index(val) for val in ptm_mapping.values()))
    cdr_ascii_index = max((ascii_letters.index(val) for val in cdr_mapping.values()))

    # jsons_sdf.head()

    df_len = len(jsons_sdf.index)
    result = {}
    for i in tqdm(range(df_len)):
        currentJson = json.loads(jsons_sdf['json_data'].iloc[i])
        for chain_type, chain_value in currentJson['ptm_predictions'].items():
            for ptm_type, ptm_value_list in chain_value.items():
                ptm_mini_name = key_get_or_create(ptm_mapping, f'{chain_type} {ptm_type}',
                                                  ascii_letters[ptm_ascii_index])
                column_ensure(dst_sdf, ptm_mini_name)
                if ptm_mini_name == ascii_letters[ptm_ascii_index]:
                    ptm_ascii_index += 1
                current_ptm = {}
                for cdr_type, cdr_ranges in currentJson['cdr_ranges'].items():
                    if f'CDR{chain_type}' in cdr_type:
                        cdr_mini_name = key_get_or_create(cdr_mapping, cdr_type, ascii_letters[cdr_ascii_index])
                        if cdr_mini_name == ascii_letters[cdr_ascii_index]:
                            cdr_ascii_index += 1
                        ptm_in_cdr = ptm_in_cdr_ranges(cdr_ranges, ptm_value_list)
                        if ptm_in_cdr:
                            current_max = max(ptm_in_cdr)
                            key_get_or_create(current_ptm, cdr_mini_name, current_max)
                            if key_get_or_create(current_ptm, 'max', current_max) < current_max:
                                current_ptm['max'] = current_max
                if current_ptm:
                    result_string = ''
                    for key, value in current_ptm.items():
                        result_string += f'{key}:{value};'

                    v_id: str = jsons_sdf.iloc[i, jsons_sdf.columns.get_loc('v_id')];
                    # jsons_sdf.iloc[i, jsons_sdf.columns.get_loc(ptm_mini_name)] = result_string[:-1]
                    v_id_dst_idx = (dst_sdf[dst_sdf['v_id'] == v_id].index.to_list()[:1] or [None])[0]
                    if v_id_dst_idx is None:
                        dst_sdf.loc[len(dst_sdf), 'v_id'] = v_id
                    dst_sdf.loc[dst_sdf['v_id'] == v_id, ptm_mini_name] = result_string[:-1]

    # jsons_sdf.head()

    # jsons_sdf.sort_values('v_id', inplace=True, ignore_index=True)
    # jsons_sdf.reindex(sorted(jsons_sdf.columns), axis=1)
    # jsons_sdf.drop(['json_data'], axis=1, inplace=True)
    # jsons_sdf.to_csv('ptm_in_cdr3_own.csv', index=False)
    #
    # with open('ptm_map3_own.json', 'wb') as f:
    #     f.write(json.dumps(ptm_mapping))
    #
    # with open('cdr_map3_own.json', 'wb') as f:
    #     f.write(json.dumps(cdr_mapping))
    #
    # df = pd.read_csv('ptm_in_cdr.csv')
    #
    # df.drop('Unnamed: 0', axis=1, inplace=True)
    #
    # df.to_feather('ptm_in_cdr.feather')
    #
    # pd.show_versions()
    pass


def build_schema(ptm_map) -> SimpleNamespace:
    # ptm_predicted = sa.Table(
    #     'ptm_predicted', meta_v2,
    #     sa.Column('v_id', sa.String(), primary_key=True, nullable=False),
    #     sa.Column('a', sa.LargeBinary, nullable=True),  # H Hydroxylysine
    #     sa.Column('b', sa.LargeBinary, nullable=True),  # H Hydroxyproline
    #     sa.Column('c', sa.LargeBinary, nullable=True),  # H Methylarginine
    #     sa.Column('d', sa.LargeBinary, nullable=True),  # H Methyllysine
    #     sa.Column('e', sa.LargeBinary, nullable=True),  # H N6-acetyllysine
    #     sa.Column('f', sa.LargeBinary, nullable=True),  # H N-linked_glycosylation
    #     sa.Column('g', sa.LargeBinary, nullable=True),  # H O-linked_glycosylation
    #     sa.Column('h', sa.LargeBinary, nullable=True),  # H Phosphoserine_Phosphothreonine
    #     sa.Column('i', sa.LargeBinary, nullable=True),  # H Phosphotyrosine
    #     sa.Column('j', sa.LargeBinary, nullable=True),  # H Pyrrolidone_carboxylic_acid
    #     sa.Column('k', sa.LargeBinary, nullable=True),  # H S-palmitoyl_cysteine
    #     sa.Column('l', sa.LargeBinary, nullable=True),  # H SUMOylation
    #     sa.Column('m', sa.LargeBinary, nullable=True),  # H Ubiquitination
    #     sa.Column('n', sa.LargeBinary, nullable=True),  # H Trpoxidation(W)
    #     sa.Column('o', sa.LargeBinary, nullable=True),  # H Asndeamidation(NGNSNT)
    #     sa.Column('p', sa.LargeBinary, nullable=True),  # H Aspisomerisation(DGDSDTDDDH)
    #     sa.Column('q', sa.LargeBinary, nullable=True),  # L Hydroxylysine
    #     sa.Column('r', sa.LargeBinary, nullable=True),  # L Hydroxyproline
    #     sa.Column('s', sa.LargeBinary, nullable=True),  # L Methylarginine
    #     sa.Column('t', sa.LargeBinary, nullable=True),  # L Methyllysine
    #     sa.Column('u', sa.LargeBinary, nullable=True),  # L N6-acetyllysine
    #     sa.Column('v', sa.LargeBinary, nullable=True),  # L N-linked_glycosylation
    #     sa.Column('w', sa.LargeBinary, nullable=True),  # L O-linked_glycosylation
    #     sa.Column('x', sa.LargeBinary, nullable=True),  # L Phosphoserine_Phosphothreonine
    #     sa.Column('y', sa.LargeBinary, nullable=True),  # L Phosphotyrosine
    #     sa.Column('z', sa.LargeBinary, nullable=True),  # L Pyrrolidone_carboxylic_acid
    #     sa.Column('A', sa.LargeBinary, nullable=True),  # L S-palmitoyl_cysteine
    #     sa.Column('B', sa.LargeBinary, nullable=True),  # L SUMOylation
    #     sa.Column('C', sa.LargeBinary, nullable=True),  # L Ubiquitination
    #     sa.Column('D', sa.LargeBinary, nullable=True),  # L Metoxidation(M)
    #     sa.Column('E', sa.LargeBinary, nullable=True),  # L Trpoxidation(W)
    #     sa.Column('F', sa.LargeBinary, nullable=True),  # L Asndeamidation(NGNSNT)
    #     sa.Column('G', sa.LargeBinary, nullable=True),  # L Aspisomerisation(DGDSDTDDDH)
    #     sa.Column('H', sa.LargeBinary, nullable=True),  # H N-linkedglycosylation(NXS/TXnotP)
    #     sa.Column('I', sa.LargeBinary, nullable=True),  # H Metoxidation(M)
    #     sa.Column('J', sa.LargeBinary, nullable=True),  # H N-terminalglutamate(VHandVL)(E)
    #     sa.Column('K', sa.LargeBinary, nullable=True),  # L N-terminalglutamate(VHandVL)(E)
    #     sa.Column('L', sa.LargeBinary, nullable=True),  # H UnpairedCys(C)
    #     sa.Column('M', sa.LargeBinary, nullable=True),  # L CD11c/CD18binding(GPR)
    #     sa.Column('N', sa.LargeBinary, nullable=True),  # L UnpairedCys(C)
    #     sa.Column('O', sa.LargeBinary, nullable=True),  # L N-linkedglycosylation(NXS/TXnotP)
    #     sa.Column('P', sa.LargeBinary, nullable=True),  # L LysineGlycation(KEKDEKED)
    #     sa.Column('Q', sa.LargeBinary, nullable=True),  # H LysineGlycation(KEKDEKED)
    #     sa.Column('R', sa.LargeBinary, nullable=True),  # L Integrinbinding(RGDRYDLDV)
    #     sa.Column('S', sa.LargeBinary, nullable=True),  # H Fragmentation(DP)
    #     sa.Column('T', sa.LargeBinary, nullable=True),  # H CD11c/CD18binding(GPR)
    #     sa.Column('U', sa.LargeBinary, nullable=True),  # L Fragmentation(DP)
    #     sa.Column('V', sa.LargeBinary, nullable=True),  # H Integrinbinding(RGDRYDLDV)
    # )

    ptm_predicted_columns = [
        sa.Column('v_id', sa.String(), primary_key=True, nullable=False),
    ]
    for (key, value) in ptm_map.items():
        col_name: str = value
        col_comment: str = key
        ptm_predicted_columns.append(sa.Column(col_name, sa.LargeBinary, nullable=True, comment=col_comment))
    ptm_predicted = sa.Table('ptm_predicted', meta_mlb, *ptm_predicted_columns)

    return SimpleNamespace(
        ptm_predicted=ptm_predicted,
    )


def load_data_ptm_predicted(conn: sae.Connection, ptm_map: dict[str, str], cdr_map: dict[str, str],
                            ptm_predicted: sa.Table):
    jsons_sdf: pd.DataFrame = pd.read_sql_table(table_name='json_files', schema=ptm_predicted.schema, con=conn);
    ptm_predicted_sdf: pd.DataFrame = pd.read_sql_table(table_name=ptm_predicted.name, schema=ptm_predicted.schema, con=conn);

    vladimir_script(jsons_sdf, ptm_map, cdr_map, ptm_predicted_sdf)

    ptm_predicted_sdf.to_sql(name=ptm_predicted.name, schema=ptm_predicted.schema, con=conn, index=False, if_exists='append')


@click.group(cls=DefaultGroup, default='main')
def cli():
    pass


@cli.command()
@click.pass_context
@click.option('--conn-str', 'conn_str',
              help="Database connection string .",
              type=str)
@click.option('--ptm-map', 'ptm_map_f',
              help='JSON file with mapping ptm names to short column names',
              type=click.File('r'))
@click.option('--cdr-map', 'cdr_map_f',
              help='JSON file with mapping ptm names to short names',
              type=click.File('r'))
def main(ctx, conn_str, ptm_map_f, cdr_map_f):
    ptm_map = json.loads(ptm_map_f.read())
    cdr_map = json.loads(cdr_map_f.read())
    s = build_schema(ptm_map)

    engine = sa.create_engine(conn_str)
    meta_mlb.create_all(bind=engine)
    with engine.begin() as conn:
        load_data_ptm_predicted(conn, ptm_map, cdr_map, s.ptm_predicted)


if __name__ == '__main__':
    cli()
