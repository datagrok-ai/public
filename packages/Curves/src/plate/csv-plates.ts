/* eslint-disable max-len */
// In a new file src/plate/csv-plates.ts, or package.ts
import * as DG from 'datagrok-api/dg';
import {Plate} from './plate';

/**
 * Parses a CSV string in a 'long' or 'tidy' format into a Plate object.
 * The CSV is expected to have a column for well positions (e.g., 'A1', 'B2').
 * @param csvContent The string content of the CSV file.
 * @param name Optional name for the plate.
 * @returns A Promise that resolves to a Plate object.
 */
export async function parsePlateFromCsv(csvContent: string, name?: string): Promise<Plate> {
  try {
    const df = DG.DataFrame.fromCsv(csvContent);
    const plate = Plate.fromTableByRow(df);

    if (name)
      plate.data.name = name;


    return plate;
  } catch (error) {
    console.error('Failed to parse CSV as a plate:', error);
    throw new Error('Could not parse the CSV file. Ensure it is in a tidy format with a well position column (e.g., "A1").');
  }
}
