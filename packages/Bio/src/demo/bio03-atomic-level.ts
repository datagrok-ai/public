import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


import {_package, toAtomicLevel} from '../package';
import $ from 'cash-dom';
import {handleError} from './utils';
import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';

const dataFn: string = 'sample/sample_FASTA.csv';

export async function demoBio03UI(): Promise<void> {
  let df: DG.DataFrame;
  let view: DG.TableView;

  try {
    await new DemoScript(
      'Atomic Level',
      'Atomic level structure of Macromolecules'
    )
      .step(`Loading Macromolecules notation 'Helm'`, async () => {
        df = await _package.files.readCsv(dataFn);
        view = grok.shell.addTableView(df);
        for (let colI: number = 0; colI < view.grid.columns.length; colI++) {
          const gCol: DG.GridColumn = view.grid.columns.byIndex(colI)!;
          if (!(['Sequence', 'Activity'].includes(gCol.name))) gCol.visible = false;
        }
      }, {
        description: `Load dataset with macromolecules of 'fasta' notation, 'PT' alphabet (protein, aminoacids).`,
        delay: 1600,
      })
      .step('To atomic level', async () => {
        const seqCol = df.getCol('Sequence');
        await toAtomicLevel(df, seqCol);
      }, {
        description: 'Get atomic level structures of Macromolecules.',
        delay: 1600,
      })
      .start();
  } catch (err: any) {
    handleError(err);
  }
}
