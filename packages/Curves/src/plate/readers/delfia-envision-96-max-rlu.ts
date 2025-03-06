import {Plate} from "../plate";
import {IPlateReader} from "../plate-reader";

export class DelfiaEnvision96MaxRluPlateReader implements IPlateReader {
  isApplicable(s: string): boolean {
    return s.startsWith('Plate information') && s.includes('DELFIA') && s.includes('Number of the wells in the plate,,,,96');
  }

  read(s: string): Plate {
    const start = s.indexOf(',01,02,03,04,05,06,07,08,09,10,11,12,');
    let end = s.indexOf('\n\n', start + 1);
    end = end != -1 ? end : s.indexOf('\r\n\r\n', start + 1);
    return Plate.fromCsv(s.substring(start, end));
  }
}