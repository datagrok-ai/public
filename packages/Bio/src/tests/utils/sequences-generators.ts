import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ALIGNMENT, ALPHABET, NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio';


export function generateManySequences(): DG.Column[] {
  const columns: DG.Column[] = [];
  columns.push(DG.Column.fromList('string', 'MSA', new Array(10 ** 6).fill('meI/hHis/Aca/N/T/dE/Thr_PO3H2/Aca/D-Tyr_Et/Tyr_ab-dehydroMe/dV/E/N/D-Orn/D-aThr//Phe_4Me')));
  columns.push(DG.Column.fromList('string', 'Activity', new Array(10 ** 6).fill('5.30751')));
  return columns;
}

export function generateLongSequence(): DG.Column[] {
  const columns: DG.Column[] = [];
  const longSequence = `meI/hHis/Aca/N/T/dE/Thr_PO3H2/Aca/D-Tyr_Et/Tyr_ab-dehydroMe/dV/E/N/D-Orn/D-aThr`.repeat(10 ** 5);
  columns.push(DG.Column.fromList('string', 'MSA', new Array(10 ** 2).fill(longSequence)));
  columns.push(DG.Column.fromList('string', 'Activity', new Array(10 ** 2).fill('7.30751')));
  return columns;
}

export function setTagsMacromolecule(col: DG.Column) {
  col.semType = DG.SEMTYPE.MACROMOLECULE;
  col.setTag(DG.TAGS.UNITS, NOTATION.SEPARATOR);
  col.setTag(bioTAGS.aligned, ALIGNMENT.SEQ_MSA);
  col.setTag(bioTAGS.alphabet, ALPHABET.UN);
  col.setTag(bioTAGS.separator, '/');
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
