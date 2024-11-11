import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package, toAtomicLevel} from '../package';
import {handleError} from './utils';
import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';
import {delay} from '@datagrok-libraries/utils/src/test';

export async function demoBio03UI(): Promise<void> {
  const dataFn: string = 'samples/HELM.csv';
  const seqColName = 'HELM';

  let df: DG.DataFrame;
  let view: DG.TableView;
  let dlg: DG.Dialog;

  try {
    await new DemoScript('Atomic Level', 'Atomic level structure of Macromolecules', false, {autoStartFirstStep: true})
      .step(`Loading Macromolecules notation 'Helm'`, async () => {
        grok.shell.windows.showContextPanel = false;
        grok.shell.windows.showProperties = false;

        df = await _package.files.readCsv(dataFn);
        view = grok.shell.addTableView(df);
        for (let colI: number = 0; colI < view.grid.columns.length; colI++) {
          const gCol: DG.GridColumn = view.grid.columns.byIndex(colI)!;
          if (!([seqColName, 'Activity'].includes(gCol.name))) gCol.visible = false;
        }
      }, {
        description: `Load dataset with macromolecules of 'fasta' notation, 'PT' alphabet (protein, aminoacids).`,
        delay: 2000,
      })
      .step('To atomic level', async () => {
        const seqCol = df.getCol(seqColName);
        await toAtomicLevel(df, seqCol, false, false);
      }, {
        description: 'Get atomic level structures of Macromolecules.',
        delay: 2000,
      })
      .step('Sketcher', async () => {
        const molColName: string = `molfile(${seqColName})`;
        df.currentCell = df.cell(1, molColName);
        const mol: string = df.currentCell.value;

        const sketcher = new DG.chem.Sketcher(DG.chem.SKETCHER_MODE.INPLACE);
        sketcher.setMolFile(mol);

        dlg = ui.dialog()
          .add(sketcher)
          .show();
        await delay(3000);
        dlg.close();
      }, {
        description: 'Display atomic level structure within a sketcher.',
        delay: 2000,
      })
      .start();
  } catch (err: any) {
    handleError(err);
  }
}
