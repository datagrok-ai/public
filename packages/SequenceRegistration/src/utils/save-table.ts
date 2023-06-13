import * as DG from 'datagrok-api/dg';
import {COL_NAMES, GENERATED_COL_NAMES, SEQUENCE_TYPES} from './const';
import {download} from './helpers';
import * as grok from 'datagrok-api/grok';
import {RegistrationSequenceParser} from './sequence-parser';

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
      if (typeCol.get(i) === SEQUENCE_TYPES.SENSE_STRAND) {
        const molfile = await getMolfile(sequenceCol.get(i), false);
        rowStr += `${molfile}\n> <Sequence>\nSense Strand\n\n`;
      } else if (typeCol.get(i) === SEQUENCE_TYPES.ANTISENSE_STRAND) {
        const molfile = await getMolfile(sequenceCol.get(i), true);
        rowStr += `${molfile}\n> <Sequence>\nAnti Sense\n\n`;
      } else if (typeCol.get(i) === SEQUENCE_TYPES.DUPLEX) {
        const obj = parser.getDuplexStrands(sequenceCol.get(i));
        const asMolfile = await getMolfile(obj.as, true);
        const as = `${asMolfile}\n> <Sequence>\nAnti Sense\n\n`;
        const ssMolfile = await getMolfile(obj.ss, false);
        const ss = `${ssMolfile}\n> <Sequence>\nSense Strand\n\n`;
        rowStr += `${await linkStrandMolfiles({senseStrands: [ss], antiStrands: [as]})}\n\n`;
      } else if ([SEQUENCE_TYPES.TRIPLEX, SEQUENCE_TYPES.DIMER].includes(typeCol.get(i))) {
        const obj = parser.getDimerStrands(sequenceCol.get(i));
        const as1Molfile = await getMolfile(obj.as1, true);
        const as1 = `${as1Molfile}\n> <Sequence>\nAnti Sense\n\n`;

        const as2Molfile = await getMolfile(obj.as2, true);
        const as2 = `${as2Molfile}\n> <Sequence>\nAnti Sense\n\n`;

        const ssMolfile = await getMolfile(obj.ss, false);
        const ss = `${ssMolfile}\n> <Sequence>\nSense Strand\n\n`;

        rowStr += `${await linkStrandMolfiles({senseStrands: [ss], antiStrands: [as1, as2]})}\n\n`;
      }

      for (const col of table.columns) {
        if (col.name !== COL_NAMES.SEQUENCE)
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

function differenceOfTwoArrays(a: any[], b: any[]): any[] {
  return a.filter((x) => !b.includes(x));
}

async function getMolfile(sequence: string, invert: boolean): Promise<string> {
  return grok.functions.call('SequenceTranslator:getMolfileFromGcrsSequence', {sequence: sequence, invert: invert});
}

async function linkStrandMolfiles(strands: { senseStrands: string[], antiStrands: string[] }) {
  return await grok.functions.call('SequenceTranslator:linkStrands', {strands: strands});
}
