import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

//import {_package} from './package-test-int';
import {before, after, category/*, expect*/, test} from '@datagrok-libraries/utils/src/test';

category('pdbGridCellRenderer', () => {
  let dfList: DG.DataFrame[];
  let viewList: DG.View[];
  let currentView: DG.ViewBase;

  before(async () => {
    dfList = [];
    viewList = [];
    currentView = grok.shell.v;
  });

  after(async () => {
    // for (const v of viewList) v.close();
    // for (const df of dfList) grok.shell.closeTable(df);
    grok.shell.v = currentView;
  });

  test('pdbDataCsv', async () => {
    const pdbDataDf: DG.DataFrame =
      await grok.dapi.files.readCsv('System:AppData/BiostructureViewer/pdb_data.csv');
    const pdbDataView: DG.TableView = grok.shell.addTableView(pdbDataDf);
    dfList.push(pdbDataDf);
    viewList.push(pdbDataView);
  });
});
