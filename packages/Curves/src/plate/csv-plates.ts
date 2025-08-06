/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import {Plate} from './plate';

// --- This is the new, powerful function used by the "Create Plate" view ---
/**
 * Parses a 'tidy' DataFrame into one or more Plate objects based on specified identifier columns.
 * @param sourceDf The source DataFrame.
 * @param plateIdColName The column that identifies the main plate group (e.g., 'barcode').
 * @param replicateIdColName An optional column that identifies technical replicates within a plate group.
 * @param plateNumberColName An optional column that identifies a plate number within a main group.
 * @returns An array of Plate objects.
 */
export function parsePlates(
  sourceDf: DG.DataFrame,
  plateIdColName: string | null,
  replicateIdColName: string | null,
  plateNumberColName: string | null
): Plate[] {
  const plates: Plate[] = [];
  const plateIdCol = plateIdColName ? sourceDf.col(plateIdColName) : null;
  const replicateIdCol = replicateIdColName ? sourceDf.col(replicateIdColName) : null;
  const plateNumberCol = plateNumberColName ? sourceDf.col(plateNumberColName) : null;

  // If no primary identifier is provided, treat the entire file as a single plate.
  if (!plateIdCol) {
    plates.push(Plate.fromTableByRow(sourceDf));
    return plates;
  }

  // Create a temporary combined key column to uniquely identify each plate/replicate/number combination.
  const comboKeyColName = '~comboKey';
  const comboKeyCol = DG.Column.string(comboKeyColName, sourceDf.rowCount).init((i) => {
    const parts = [plateIdCol.get(i)];
    if (plateNumberCol) parts.push(plateNumberCol.get(i));
    if (replicateIdCol) parts.push(replicateIdCol.get(i));
    return parts.join('_');
  });

  const uniqueIds = comboKeyCol.categories;

  if (uniqueIds.length > 1) {
    // If multiple unique combinations are found, split the DataFrame into separate plates.
    for (const id of uniqueIds) {
      const rowMask = DG.BitSet.create(sourceDf.rowCount, (i) => comboKeyCol.get(i) === id);
      const plateDf = sourceDf.clone(rowMask);

      const plate = Plate.fromTableByRow(plateDf);
      plate.barcode = id; // The plate's barcode becomes the combined unique ID.
      plates.push(plate);
    }
    return plates;
  }

  // Fallback for a single plate (or a file with only one unique combo-key).
  const plate = Plate.fromTableByRow(sourceDf);
  if (uniqueIds.length === 1)
    plate.barcode = uniqueIds[0];

  plates.push(plate);
  return plates;
}


// --- This is the backward-compatibility wrapper for other parts of the app ---
/**
 * Parses a CSV string using default settings for simple file handlers.
 * It defaults to looking for a 'barcode' column to split by, preserving the original behavior.
 * @param csvContent The string content of the CSV file.
 * @returns A Promise that resolves to an array of Plate objects.
 */
export async function parsePlateFromCsv(csvContent: string): Promise<Plate[]> {
  try {
    const df = DG.DataFrame.fromCsv(csvContent);
    // Call our new, powerful function with default parameters ('barcode' for plate, null for others).
    return parsePlates(df, 'barcode', null, null);
  } catch (error) {
    console.error('Failed to parse CSV as a plate:', error);
    throw new Error('Could not parse the CSV file. Ensure it is in a tidy format with a well position column (e.g., "A1").');
  }
}
