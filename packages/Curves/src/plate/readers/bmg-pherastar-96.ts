import {Plate} from "../plate";
import {IPlateReader} from "../plate-reader";
import { parseExcelPosition } from "../utils";

export class BmgPherastar96PlateReader implements IPlateReader {
  isApplicable(s: string): boolean {
    return s.includes('No. of Channels / Multichromatics:') && s.includes('No. of Cycles:');
  }

  read(s: string): Plate {
    // Find the data section starting with "A01;"
    const dataStart = s.indexOf('A01;');
    if (dataStart === -1) throw new Error('Data section not found');

    // Split into lines and process each line
    const lines = s.substring(dataStart).split('\n')
      .map(line => line.trim())
      .filter(line => line.length > 0)
      .map(line => line.split(';').map(s => s.trim()));

    const plate = Plate.autoSize(lines.map((line) => parseExcelPosition(line[0])));

    // Create columns for both channels
    const channelA = plate.data.columns.addNewFloat('Channel A');
    const channelB = plate.data.columns.addNewFloat('Channel B');

    // Process each line
    for (const line of lines) {
      const [pos, valA, valB] = line;
      if (!pos || pos.length !== 3) continue;

      const [row, col] = parseExcelPosition(pos);

      // Set values in the plate
      const idx = plate._idx(row, col);
      channelA.set(idx, parseFloat(valA), false);
      channelB.set(idx, parseFloat(valB), false);
    }

    return plate;
  }
}