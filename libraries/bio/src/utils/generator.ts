import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {ALPHABET, NOTATION, TAGS} from './macromolecule';

export function generateManySequences(): DG.Column[] {
  const columns: DG.Column[] = [];
  columns.push(DG.Column.fromList('string', 'MSA',
    new Array(10 ** 6).fill(
      'meI/hHis/Aca/N/T/dE/Thr_PO3H2/Aca/D-Tyr_Et/Tyr_ab-dehydroMe/dV/E/N/D-Orn/D-aThr//Phe_4Me')),
  );
  columns.push(DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Activity', new Array(10 ** 6).fill(5.30751)));
  return columns;
}

/** Generates the column 'MSA' with sequences length of order 10^6 and the 'Activity' float column. */
export function generateLongSequence(length: number = 10 ** 5): DG.Column[] {
  const columns: DG.Column[] = [];
  const longSequence =
    `meI/hHis/Aca/N//dE/Thr_PO3H2/Aca/D-Tyr_Et/Tyr_ab-dehydroMe`.repeat(Math.ceil(length / 10)).slice(0, -1);
  const msaCol = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'MSA', new Array(10 ** 2).fill(longSequence));
  msaCol.semType = DG.SEMTYPE.MACROMOLECULE;
  msaCol.setTag(DG.TAGS.UNITS, NOTATION.SEPARATOR);
  msaCol.setTag(TAGS.separator, '/');
  msaCol.setTag(TAGS.alphabet, ALPHABET.UN);
  msaCol.setTag(TAGS.alphabetIsMultichar, 'true');

  columns.push(msaCol);
  columns.push(DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Activity', new Array(10 ** 2).fill(7.30751)));
  return columns;
}
