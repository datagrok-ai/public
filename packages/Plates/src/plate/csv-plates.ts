/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import {Plate} from './plate';

export interface PlateFile {
 plate: Plate;
 commonProperties: Map<string, any>;
 reconciliationMap: Map<string, string>;
}

function findCommonProperties(plateDf: DG.DataFrame): Map<string, any> {
  const commonProperties = new Map<string, any>();
  const IGNORED_PROPERTIES = new Set<string>([
  ]);

  for (const col of plateDf.columns) {
    if (!IGNORED_PROPERTIES.has(col.name) && col.stats.uniqueCount === 1) {
      const value = col.get(0);
      if (value !== null && value !== '')
        commonProperties.set(col.name, value);
    }
  }
  return commonProperties;
}

export function parsePlates(
  sourceDf: DG.DataFrame,
  identifierColNames: (string | null)[]
): Array<{ plate: Plate, commonProperties: Map<string, any> }> {
  const results: { plate: Plate, commonProperties: Map<string, any> }[] = [];

  const validIdentifierNames = identifierColNames.filter((name): name is string =>
    name !== null && sourceDf.columns.contains(name));

  if (validIdentifierNames.length === 0) {
    const plate = Plate.fromTableByRow(sourceDf);
    const commonProperties = findCommonProperties(sourceDf);
    results.push({plate, commonProperties});
    return results;
  }

  const identifierCols = validIdentifierNames.map((name) => sourceDf.col(name)!);
  const comboKeyColName = '~comboKey';
  const comboKeyCol = DG.Column.string(comboKeyColName, sourceDf.rowCount).init((i) => {
    return identifierCols.map((col) => col.get(i)).join('_');
  });

  const uniqueIds = comboKeyCol.categories;

  for (const id of uniqueIds) {
    const rowMask = DG.BitSet.create(sourceDf.rowCount, (i) => comboKeyCol.get(i) === id);
    const plateDf = sourceDf.clone(rowMask);

    const plate = Plate.fromTableByRow(plateDf);
    plate.barcode = id;
    const commonProperties = findCommonProperties(plateDf);
    results.push({plate, commonProperties});
  }

  return results;
}


export async function parsePlateFromCsv(csvContent: string): Promise<Plate[]> {
  try {
    const df = DG.DataFrame.fromCsv(csvContent);
    return parsePlates(df, ['Destination Plate Barcode']).map((res) => res.plate);
  } catch (error) {
    console.error('Failed to parse CSV as a plate:', error);
    throw new Error('Could not parse the CSV file. Ensure it is in a tidy format with a well position column (e.g., "A1").');
  }
}
