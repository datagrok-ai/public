import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

export async function searchSubstructure(molStr: string, dbdf: DG.DataFrame): Promise<DG.DataFrame | null> {
  const bitset = await grok.chem.searchSubstructure(dbdf.getCol('molecule'), molStr);

  if (bitset === null)
    return null;

  dbdf.filter.copyFrom(bitset);
  return dbdf;
}

export async function findSimilar(molStr: string, limit: number, cutoff: number, dbdf: DG.DataFrame,
): Promise<DG.DataFrame | null> {
  const searchdf = await grok.chem.findSimilar(dbdf.getCol('molecule'), molStr, {'limit': limit, 'cutoff': cutoff});

  if (searchdf == null)
    return null;

  if (dbdf.col('index') === null)
    dbdf.columns.addNewInt('index').init((i) => i);

  return dbdf.join(searchdf, ['index'], ['index'], ['molecule'], ['molecule'], DG.JOIN_TYPE.INNER, true);
}
