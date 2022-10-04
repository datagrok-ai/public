import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {DataFrame} from 'datagrok-api/dg';

import {NOTATION, ALIGNMENT, UnitsHandler, ALPHABET} from '@datagrok-libraries/bio/src/utils/units-handler'

export function generateManySequences(): DG.Column[] {
  let columns: DG.Column[] = [];
  columns.push(DG.Column.fromList('string', 'MSA', new Array(10 ** 6).fill('meI/hHis/Aca/N/T/dE/Thr_PO3H2/Aca/D-Tyr_Et/Tyr_ab-dehydroMe/dV/E/N/D-Orn/D-aThr//Phe_4Me')));
  columns.push(DG.Column.fromList('string', 'Activity', new Array(10 ** 6).fill('5.30751')));
  return columns;
}

export function generateLongSequence(): DG.Column[] {
  let columns: DG.Column[] = [];
  const longSequence = `meI/hHis/Aca/N/T/dE/Thr_PO3H2/Aca/D-Tyr_Et/Tyr_ab-dehydroMe/dV/E/N/D-Orn/D-aThr`.repeat(10 ** 5);
  columns.push(DG.Column.fromList('string', 'MSA', new Array(10 ** 2).fill(longSequence)));
  columns.push(DG.Column.fromList('string', 'Activity', new Array(10 ** 2).fill('7.30751')));
  return columns;
}

export function setTagsMacromolecule(col: DG.Column) {
  col.semType = DG.SEMTYPE.MACROMOLECULE;
  col.setTag(DG.TAGS.UNITS, NOTATION.SEPARATOR);
  col.setTag(UnitsHandler.TAGS.aligned, ALIGNMENT.SEQ_MSA);
  col.setTag(UnitsHandler.TAGS.alphabet, ALPHABET.UN);
  col.setTag(UnitsHandler.TAGS.separator, '/');
  return col;
}

export function performanceTest(generateFunc: () => DG.Column[], testName: string) {
  const columns = generateFunc();
  const df: DG.DataFrame = DG.DataFrame.fromColumns(columns);
  const startTime: number = Date.now();
  const col: DG.Column = df.columns.byName('MSA');
  setTagsMacromolecule(col);
  grok.shell.addTableView(df);

  const endTime: number = Date.now();
  const elapsedTime: number = endTime - startTime;
  console.log(`Performance test: ${testName}: ${elapsedTime}ms`);
}
