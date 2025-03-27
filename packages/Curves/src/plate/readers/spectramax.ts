import * as DG from 'datagrok-api/dg';
import { Plate } from "../plate";
import { IPlateReader } from "../plate-reader";

export class SpectramaxPlateReader implements IPlateReader {
  isApplicable(s: string): boolean {
    return s.includes('Plate:') && s.includes('##BLOCKS=');
  }

  read(content: string): Plate {
    const lines = content.split(/\r\n|\r|\n/);
    const header = lines[2].split('\t').filter(s => s !== '');
    const plateCount = parseInt(lines[1].split('\t')[8]);

    // empty measurements names in the beginning are ok - need for indexing, will be ignored later
    const plateMeasurementNames = header.filter((s, i) => s == "" && i < 5 || !/^-?\d+$/.test(s));
    const cols = header.filter((s) => /^-?\d+$/.test(s)).length;
    const rows = cols * 2 / 3;
    const plates = [];

    for (let plateIdx = 0; plateIdx < plateCount; plateIdx++) {
      const plate = new Plate(rows, cols);
      const headerValues = lines[plateIdx * (rows + 1) + 3].split('\t');
      const plateMeasurements = plateMeasurementNames
        .map((v, i) => v + ': ' + headerValues[i])
        .filter((v, i) => v != '')
        .join(', ');

      const data = plate.data.columns.addNewFloat(plateMeasurements);
      for (let row = 0; row < rows; row++) {
        const rowIdx = plateIdx * (rows + 1) + 3 + row;
        const values = lines[rowIdx].split('\t');
        for (let col = 0; col < cols; col++) {
          data.set(row * cols + col, parseFloat(values[col + 2].replace(',', '.')), false);
        }
      }

      plates.push(plate);
    }

    return Plate.fromPlates(plates);
  }
}