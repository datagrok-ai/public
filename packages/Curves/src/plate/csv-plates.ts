/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import {Plate} from './plate';

export interface PlateFile {
  plate: Plate;
  commonProperties: Map<string, any>;
  reconciliationMap: Map<string, string>;
}

/**
 * Analyzes a DataFrame for a single plate to find columns with constant values.
 * @param plateDf A DataFrame containing data for only one plate.
 * @returns A Map of common property names to their values.
 */
function findCommonProperties(plateDf: DG.DataFrame): Map<string, any> {
  const commonProperties = new Map<string, any>();
  const IGNORED_PROPERTIES = new Set([
    'Well Position', 'Robust Inhibition (%)', 'Concentration (M)', 'Signal Emission',
    'Destination Plate Barcode', 'Plate Number', 'Technical Duplicate ID', 'Vial Barcode',
    'Automatic Data Point Exclusion (\'Excluded\', NA)', 'Manual Data Point Exclusion (\'Excluded\', NA)',
    'SMILES', 'Compound Iso ID'
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

/**
 * Parses a 'tidy' DataFrame and returns an array of objects, each containing a Plate and its common properties.
 */
export function parsePlates(
  sourceDf: DG.DataFrame,
  plateIdColName: string | null,
  replicateIdColName: string | null,
  plateNumberColName: string | null
): { plate: Plate, commonProperties: Map<string, any> }[] { // <-- Return an intermediate object
  const results: { plate: Plate, commonProperties: Map<string, any> }[] = [];
  const plateIdCol = plateIdColName ? sourceDf.col(plateIdColName) : null;

  if (!plateIdCol) {
    const plate = Plate.fromTableByRow(sourceDf);
    const commonProperties = findCommonProperties(sourceDf);
    results.push({plate, commonProperties});
    return results;
  }

  const replicateIdCol = replicateIdColName ? sourceDf.col(replicateIdColName) : null;
  const plateNumberCol = plateNumberColName ? sourceDf.col(plateNumberColName) : null;

  const comboKeyColName = '~comboKey';
  const comboKeyCol = DG.Column.string(comboKeyColName, sourceDf.rowCount).init((i) => {
    const parts = [plateIdCol.get(i)];
    if (plateNumberCol) parts.push(plateNumberCol.get(i));
    if (replicateIdCol) parts.push(replicateIdCol.get(i));
    return parts.join('_');
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


/**
 * Backward-compatibility wrapper.
 * @returns A Promise that resolves to an array of Plate objects.
 */
export async function parsePlateFromCsv(csvContent: string): Promise<Plate[]> {
  try {
    const df = DG.DataFrame.fromCsv(csvContent);
    // Note: The return type here is Plate[] for backward compatibility.
    // New code should use parsePlates directly.
    return parsePlates(df, 'Destination Plate Barcode', null, null).map((res) => res.plate);
  } catch (error) {
    console.error('Failed to parse CSV as a plate:', error);
    throw new Error('Could not parse the CSV file. Ensure it is in a tidy format with a well position column (e.g., "A1").');
  }
}
