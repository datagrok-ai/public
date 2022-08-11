import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {drugBankSearchResult} from './utils';

export async function drugBankSubstructureSearch(
  molString: string, dbdf: DG.DataFrame): Promise<drugBankSearchResult> {
  const bitset = await grok.chem.searchSubstructure(dbdf!.getCol('molecule'), molString);

  if (bitset === null)
    return null;

  dbdf!.filter.copyFrom(bitset);
  return dbdf!;
}

export async function drugBankSimilaritySearch(
  molString: string, limit: number, cutoff: number, dbdf: DG.DataFrame): Promise<drugBankSearchResult> {
  const searchdf = await grok.chem.findSimilar(dbdf!.getCol('molecule'), molString, {'limit': limit, 'cutoff': cutoff});

  if (searchdf == null)
    return null;

  if (dbdf!.col('index') === null)
    (dbdf!.columns as DG.ColumnList).addNewInt('index').init((i) => i);

  return dbdf!.join(searchdf, ['index'], ['index'], ['molecule'], ['molecule'], DG.JOIN_TYPE.INNER, true);
}
