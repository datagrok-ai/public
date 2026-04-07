import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {ALIGNMENT, ALPHABET, NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {expect} from '@datagrok-libraries/test/src/test';

import {awaitGrid} from '../utils';

export async function performanceTest(generateFunc: () => DG.Column[], testName: string) {
  const columns = generateFunc();
  const df: DG.DataFrame = DG.DataFrame.fromColumns(columns);
  await grok.data.detectSemanticTypes(df);
  const startTime: number = Date.now();
  const col: DG.Column = df.columns.byName('MSA');
  const view: DG.TableView = grok.shell.addTableView(df);

  await awaitGrid(view.grid);
  expect(view.grid.dataFrame.id, df.id);

  const endTime: number = Date.now();
  const elapsedTime: number = endTime - startTime;
  console.log(`Performance test: ${testName}: ${elapsedTime}ms`);
}
