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

export function parsePlates(
  sourceDf: DG.DataFrame,
  identifierColNames: (string | null)[]
): { plate: Plate, commonProperties: Map<string, any> }[] {
// END of CHANGE
  const results: { plate: Plate, commonProperties: Map<string, any> }[] = [];

  // START of CHANGE: New logic to handle dynamic identifiers
  const validIdentifierNames = identifierColNames.filter((name): name is string =>
    name !== null && sourceDf.columns.contains(name));

  if (validIdentifierNames.length === 0) {
    // If no valid identifiers are provided, treat the entire file as a single plate.
    const plate = Plate.fromTableByRow(sourceDf);
    const commonProperties = findCommonProperties(sourceDf);
    results.push({plate, commonProperties});
    return results;
  }

  const identifierCols = validIdentifierNames.map((name) => sourceDf.col(name)!);
  const comboKeyColName = '~comboKey';
  const comboKeyCol = DG.Column.string(comboKeyColName, sourceDf.rowCount).init((i) => {
    // Create a unique key by joining values from all specified identifier columns
    return identifierCols.map((col) => col.get(i)).join('_');
  });
  // END of CHANGE

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
    // Note: The return type here is Plate[] for backward compatibility.
    // New code should use parsePlates directly.
    // For compatibility, we hardcode the old default. This function is likely deprecated anyway.
    return parsePlates(df, ['Destination Plate Barcode']).map((res) => res.plate);
  } catch (error) {
    console.error('Failed to parse CSV as a plate:', error);
    throw new Error('Could not parse the CSV file. Ensure it is in a tidy format with a well position column (e.g., "A1").');
  }
}
