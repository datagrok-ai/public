import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {drugBankSearchResult} from './utils';

export async function drugBankSubstructureSearch(
  molString: string, substructLibrary: boolean, dbdf: DG.DataFrame): Promise<drugBankSearchResult> {
  const bitset = await grok.chem.searchSubstructure(
    dbdf!.getCol('SMILES'), molString, {'substructLibrary': substructLibrary});

  if (bitset === null)
    return null;

  dbdf!.filter.copyFrom(bitset);
  return dbdf!;
}

export async function drugBankSimilaritySearch(
  molString: string, limit: number, cutoff: number, dbdf: DG.DataFrame): Promise<drugBankSearchResult> {
  const searchdf = await grok.chem.findSimilar(dbdf!.getCol('SMILES'), molString, {'limit': limit, 'cutoff': cutoff});

  if (searchdf == null)
    return null;

  if (dbdf!.col('index') === null)
    (dbdf!.columns as DG.ColumnList).addNewInt('index').init((i) => i);

  return dbdf!.join(searchdf, ['index'], ['index'], ['SMILES'], ['SMILES'], DG.JOIN_TYPE.INNER, true);
}
