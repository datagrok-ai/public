/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Subscription} from 'rxjs';

export async function _testDetectors(path: string, detector: DG.Func): Promise<DG.DataFrame> {
  const pi = DG.TaskBarProgressIndicator.create('Test detectors...');

  const fileList = await grok.dapi.files.list(path, true, '');
  const csvList = fileList.filter((fi) => fi.fileName.endsWith('.csv'));

  
  let readyCount = 0;
  const res = [];

  for (const fileInfo of csvList) {
    try{
      const csv = await grok.dapi.files.readAsText(path + fileInfo.fullPath);
      const df = DG.DataFrame.fromCsv(csv);
      for (const col of df.columns) {
        const semType: string | null = await detector.apply({col: col});
        if (semType !== null) {
          res.push({
            file: fileInfo.path, result: 'detected', column: col.name,
            message: `semType is ${semType}`
          });
        }
      }
    }
    catch (err: unknown) {
      res.push({
        file: fileInfo.path, result: 'error', column: null,
        message: err instanceof Error ? err.message : (err as Object).toString(),
      });
    }
    finally {
      readyCount += 1;
      pi.update(100 * readyCount / csvList.length, `Test ${fileInfo.fileName}`);
    }
  }
  grok.shell.info(`Test ${path} for ${detector.name} finished.`);
  pi.close();
  const resDf = DG.DataFrame.fromObjects(res)!;
  resDf.name = `datasets_${detector.name}_${path}`;
  return resDf;
}

let testDialog: DG.Dialog | null = null;
let testDialogSubs: Subscription[] = [];

export function _testDetectorsDialog(): void {
  const funcArray = DG.Func.find({tags: ['semTypeDetector']});
  // TODO: consider automatic choice of connections
  const dirsArray = ['Demo:Files/', 'System:AppData/'];

  const pathInput = ui.choiceInput('Path', dirsArray[0], dirsArray);
  const detectorInput = ui.choiceInput('Detector', funcArray[0], funcArray);

  if (testDialog == null) {
    testDialog = ui.dialog('Test semType detectors')
      .add(ui.div([
        pathInput.root,
        detectorInput.root
      ]))
      .onOK(async () => {
        const path = pathInput.value;
        const detector: DG.Func = detectorInput.value;
        const df = await _testDetectors(path, detector);
        grok.shell.addTableView(df);
      })
      .show();

    testDialogSubs.push(testDialog.onClose.subscribe((value) => {
      testDialogSubs.forEach((s) => { s.unsubscribe(); });
      testDialogSubs = [];
      testDialog = null;
    }));
  }
}
