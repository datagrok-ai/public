import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

export function generateManySequences(): string {
  let csvData = `MSA,Activity
meI/hHis/Aca/N/T/dE/Thr_PO3H2/Aca/D-Tyr_Et/Tyr_ab-dehydroMe/dV/E/N/D-Orn/D-aThr//Phe_4Me,5.30751`;
  for (let i = 0; i < 10 ** 6; i++) {
    csvData += `\n meI/hHis/Aca/N/T/dE/Thr_PO3H2/Aca/D-Tyr_Et/Tyr_ab-dehydroMe/dV/E/N/D-Orn/D-aThr//Phe_4Me,5.30751`;
  }
  return csvData;
}

export function generateLongSequence(): string {
  let longSequence = `meI/hHis/Aca/N/T/dE/Thr_PO3H2/Aca/D-Tyr_Et/Tyr_ab-dehydroMe/dV/E/N/D-Orn/D-aThr`;
  for (let i = 0; i < 10 ** 5; i++) {
    longSequence += `/Aca/N/T/dE/Thr_PO3H2/Aca/D-Tyr_Et/Tyr_ab-dehydroMe/dV/dv`;
  }
  longSequence += `//Phe_4Me,5.30751`;
  let csvData = `MSA,Activity `;
  for (let i = 0; i <= 10 ** 1 * 4; i++) {
    csvData += `\n ${longSequence}`;
  }
  return csvData;
}
export function setTagsMacromolecule(col: DG.Column) {
  col.semType = DG.SEMTYPE.MACROMOLECULE;
  col.setTag('units', 'separator');
  col.setTag('aligned', 'SEQ.MSA');
  col.setTag('alphabet', 'UN');
  col.setTag('separator', '/');
  return col;
}

export function performanceTest(generateFunc: () => string,testName: string) {
  const startTime: number = Date.now();
  const csv = generateFunc();
  const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
  const col: DG.Column = df.columns.byName('MSA');
  setTagsMacromolecule(col);
  grok.shell.addTableView(df);

  const endTime: number = Date.now();
  const elapsedTime: number = endTime - startTime;
  console.log(`Performance test: ${testName}: ${elapsedTime}ms`);
}
