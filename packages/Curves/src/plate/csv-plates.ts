// In src/plate/csv-plates.ts

/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import {Plate} from './plate';

/**
 * Parses a CSV string in a 'long' or 'tidy' format into one or more Plate objects.
 * If a 'barcode' column with multiple unique values is found, it splits the data
 * into a separate plate for each barcode.
 * @param csvContent The string content of the CSV file.
 * @returns A Promise that resolves to an array of Plate objects.
 */
export async function parsePlateFromCsv(csvContent: string): Promise<Plate[]> {
  try {
    const df = DG.DataFrame.fromCsv(csvContent);
    const plates: Plate[] = [];

    const barcodeColName = df.columns.names().find((name) => name.toLowerCase() === 'barcode');
    if (barcodeColName) {
      const barcodeCol = df.col(barcodeColName)!;
      const uniqueBarcodes = barcodeCol.categories;

      if (uniqueBarcodes.length > 1) {
        for (const barcode of uniqueBarcodes) {
          const rowMask = DG.BitSet.create(df.rowCount, (i) => barcodeCol.get(i) === barcode);
          const plateDf = df.clone(rowMask);

          const plate = Plate.fromTableByRow(plateDf);
          plate.barcode = barcode;
          plates.push(plate);
        }
        return plates;
      }
    }
    const plate = Plate.fromTableByRow(df);
    plates.push(plate);
    return plates;
  } catch (error) {
    console.error('Failed to parse CSV as a plate:', error);
    throw new Error('Could not parse the CSV file. Ensure it is in a tidy format with a well position column (e.g., "A1").');
  }
}
