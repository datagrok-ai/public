/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {NotationConverter} from '@datagrok-libraries/bio/src/utils/notation-converter';
import {_package} from '../package';

export async function _runEnumerator(): Promise<void> {
  const df: DG.DataFrame = await _package.files.readCsv('samples/enumerator.csv');
  const separatorCol = df.col('Synthesis');
  separatorCol?.setTag('quality', 'Macromolecule');
  separatorCol?.setTag('units', 'separator');
  separatorCol?.setTag('separator', '-');
  separatorCol?.setTag('cell.renderer', 'sequence');
  separatorCol?.setTag('aligned', 'SEQ');
  separatorCol?.setTag('alphabet', 'PT');
  const tableView = DG.TableView.create(df);
  grok.shell.addView(tableView);
}

async function enumerator(molColumn: DG.Column): Promise<void> {
  function getCyclicHelm(inputHelm: string): string {
    const leftWrapper = 'PEPTIDE1{';
    const rightWrapper = '}$$$$';
    const seq = inputHelm.replace(leftWrapper, '').replace(rightWrapper, '');
    const lastMonomerNumber = seq.split('.').length;
    const result = inputHelm.replace(rightWrapper,
      `}$PEPTIDE1,PEPTIDE1,${lastMonomerNumber}:R2-1:R1${'$'.repeat(6)}V2.0`);
    console.log('result:', result);
    return result;
  }

  const df = molColumn.dataFrame;
  const nc = new NotationConverter(molColumn);
  const sourceHelmCol = nc.convert(NOTATION.HELM);
  df.columns.add(sourceHelmCol);
  const targetList = sourceHelmCol.toList().map((helm) => getCyclicHelm(helm));
  const targetHelmCol = DG.Column.fromList('string', 'Enumerator(helm)', targetList);
  targetHelmCol.setTag('quality', 'Macromolecule');
  targetHelmCol.setTag('units', 'helm');
  targetHelmCol.setTag('cell.renderer', 'helm');
  targetHelmCol.setTag('aligned', 'SEQ');
  targetHelmCol.setTag('alphabet', 'PT');
  df.columns.add(targetHelmCol);
  await grok.data.detectSemanticTypes(df);
}

export function getEnumeratorWidget(molColumn: DG.Column): DG.Widget {
  const btn = ui.bigButton('Run', async () => enumerator(molColumn));

  const div = ui.div([btn]);

  return new DG.Widget(div);
}
