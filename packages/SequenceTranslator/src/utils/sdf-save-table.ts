import * as DG from 'datagrok-api/dg';
import {COL_NAMES, GENERATED_COL_NAMES, SEQUENCE_TYPES} from '../autostart/constants';
import {differenceOfTwoArrays, download} from '../utils/helpers';
import * as grok from 'datagrok-api/grok';
import {SYNTHESIZERS} from '../model/const';
import {sequenceToMolV3000} from '../utils/structures-works/from-monomers';
import {RegistrationSequenceParser} from '../model/registration-sequence-parser';
import {linkStrandsV3000} from '../utils/structures-works/mol-transformations';

export async function sdfSaveTable(table: DG.DataFrame, onError: (rowI: number, err: any) => void) {
  if (GENERATED_COL_NAMES.some((colName) => !table.columns.contains(colName))) {
    const absentColNames = differenceOfTwoArrays(GENERATED_COL_NAMES, table.columns.names()).join(`', '`);
    grok.shell.warning(`File saved without columns '${absentColNames}'`);
  }

  const sequenceCol = table.getCol(COL_NAMES.SEQUENCE);
  const typeCol = table.getCol(COL_NAMES.TYPE);

  let resultStr = '';
  const rowCount = table.rowCount;
  for (let i = 0; i < rowCount; i++) {
    try {
      let rowStr = '';
      const parser = new RegistrationSequenceParser();
      const format = SYNTHESIZERS.GCRS; //getFormat(sequenceCol.get(i))!;
      if (typeCol.get(i) == SEQUENCE_TYPES.SENSE_STRAND) {
        rowStr += `${sequenceToMolV3000(sequenceCol.get(i), false, format)}\n> <Sequence>\nSense Strand\n\n`;
      } else if (typeCol.get(i) == SEQUENCE_TYPES.ANTISENSE_STRAND) {
        rowStr += `${sequenceToMolV3000(sequenceCol.get(i), true, format)}\n> <Sequence>\nAnti Sense\n\n`;
      } else if (typeCol.get(i) == SEQUENCE_TYPES.DUPLEX) {
        const obj = parser.getDuplexStrands(sequenceCol.get(i));
        const as = `${sequenceToMolV3000(obj.as, true, format)}\n> <Sequence>\nAnti Sense\n\n`;
        const ss = `${sequenceToMolV3000(obj.ss, false, format)}\n> <Sequence>\nSense Strand\n\n`;
        rowStr += `${linkStrandsV3000({senseStrands: [ss], antiStrands: [as]}, true)}\n\n`;
      } else if ([SEQUENCE_TYPES.TRIPLEX, SEQUENCE_TYPES.DIMER].includes(typeCol.get(i))) {
        const obj = parser.getDimerStrands(sequenceCol.get(i));
        const as1 = `${sequenceToMolV3000(obj.as1, true, format)}\n> <Sequence>\nAnti Sense\n\n`;
        const as2 = `${sequenceToMolV3000(obj.as2, true, format)}\n> <Sequence>\nAnti Sense\n\n`;
        const ss = `${sequenceToMolV3000(obj.ss, false, format)}\n> <Sequence>\nSense Strand\n\n`;
        rowStr += `${linkStrandsV3000({senseStrands: [ss], antiStrands: [as1, as2]}, true)}\n\n`;
      }

      for (const col of table.columns) {
        if (col.name != COL_NAMES.SEQUENCE)
          rowStr += `> <${col.name}>\n${col.get(i)}\n\n`;
      }

      rowStr += '$$$$\n';

      resultStr += rowStr;
    } catch (err: any) {
      onError(i, err);
    }
  }

  download(`${table.name}.sdf`, encodeURIComponent(resultStr));
}
