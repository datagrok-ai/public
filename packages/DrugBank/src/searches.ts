import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {COLUMN_NAMES} from './widgets';

export async function searchSubstructure(molStr: string, dbdf: DG.DataFrame): Promise<DG.DataFrame | null> {
  const bitset = await grok.chem.searchSubstructure(dbdf.getCol(COLUMN_NAMES.MOLECULE), molStr);

  if (bitset === null)
    return null;

  dbdf.filter.copyFrom(bitset);
  return dbdf;
}

export async function findSimilar(molStr: string, limit: number, cutoff: number, dbdf: DG.DataFrame,
): Promise<DG.DataFrame | null> {
  const searchdf = await grok.chem.findSimilar(dbdf.getCol(COLUMN_NAMES.MOLECULE), molStr, {limit, cutoff});

  if (searchdf === null)
    return null;

  const indexes = searchdf.getCol('index').getRawData();
  const idCol: DG.Column<string> = dbdf.getCol(COLUMN_NAMES.DRUGBANK_ID);
  const resultIdCol: DG.Column<string> = searchdf.columns.addNew(idCol.name, idCol.type);
  const nameCol: DG.Column<string> = dbdf.getCol(COLUMN_NAMES.COMMON_NAME);
  const resultNameCol: DG.Column<string> = searchdf.columns.addNew(nameCol.name, nameCol.type);
  for (let idx = 0; idx < indexes.length; idx++) {
    const piv = indexes[idx];
    resultIdCol.set(idx, idCol.get(piv));
    resultNameCol.set(idx, nameCol.get(piv));
  }

  return searchdf;
}
