import * as DG from 'datagrok-api/dg';
import {Plate} from "../plate";
import {IPlateReader} from "../plate-reader";

export class DelfiaEnvisionByRowPlateReader implements IPlateReader {
  isApplicable(s: string): boolean {
    return s.includes('PltRpt=1') && s.includes('Group=1');
  }

  read(content: string): Plate {
    let end = content.indexOf('\n\n');
    end = end == -1 ? content.indexOf('\r\n\r\n') : end;
    const csv = content.substring(content.indexOf('\n') + 1, end);
    const table = DG.DataFrame.fromCsv(csv, {delimiter: '\t'});

    return Plate.fromTableByRow(table);
  }
}