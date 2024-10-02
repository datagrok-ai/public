import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import {fromEvent} from 'rxjs';

import {
  after, before, category, test, expect, awaitCheck, delay, testEvent
} from '@datagrok-libraries/utils/src/test';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {
  getUserLibSettings, setUserLibSettings
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';

import {awaitGrid, initHelmMainPackage} from './utils';
import {TAGS as helmTAGS} from '../constants';

import {_package} from '../package-test';

category('renderers', () => {
  let monomerLibHelper: IMonomerLibHelper;
  /** Backup actual user's monomer libraries settings */
  let userLibSettings: UserLibSettings;

  before(async () => {
    await initHelmMainPackage();

    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();

    // Test 'helm' requires default monomer library loaded
    await monomerLibHelper.loadMonomerLibForTests(); // load default libraries
  });

  after(async () => {
    // UserDataStorage.put() replaces existing data
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true); // load user settings libraries
  });


  test('missedInLib', async () => { await _testMissedInLib(); });

  test('scatterPlotTooltip', async () => { await _testScatterPlotTooltip(); });

  async function _testMissedInLib() {
    const df: DG.DataFrame = await grok.dapi.files.readCsv('System:AppData/Helm/tests/sample_HELM-missedInLib.csv');
    const helmCol: DG.Column<string> = df.getCol('HELM');

    const tv: DG.TableView = grok.shell.addTableView(df);
    await awaitCheck(() => {
      return $(tv.root).find('.d4-grid canvas').length > 0;
    }, 'Table view canvas not found', 100);

    expect(helmCol.semType, DG.SEMTYPE.MACROMOLECULE);
    expect(helmCol.meta.units, NOTATION.HELM);
    expect(helmCol.getTag(DG.TAGS.CELL_RENDERER), 'helm');

    const cellRendererErrorJson: string = helmCol.getTag(helmTAGS.cellRendererRenderError);
    if (!!cellRendererErrorJson) {
      const cellRendererError: any = JSON.parse(cellRendererErrorJson);
      const err = new Error(cellRendererError['message']);
      err.stack = cellRendererError['stack'];
      throw err;
    }
  }

  const helmCoordsCsv = `seq,x,y
PEPTIDE1{I.H.A.N.T.Thr_PO3H2.Aca.D-Tyr_Et}$$$$,0,0
RNA1{r(A)p.r(C)p.r(G)p.r(U)p}$$$$,1,0
RNA1{d(A)p.d(C)p.d(G)p.d(U)p}|PEPTIDE1{I.H.A.N.T.Thr_PO3H2}$$$$,0,1
`;

  async function _testScatterPlotTooltip(): Promise<void> {
    const df = DG.DataFrame.fromCsv(helmCoordsCsv);
    df.currentRowIdx = 0;
    const view = grok.shell.addTableView(df);
    const sp: DG.ScatterPlotViewer = df.plot.scatter({x: 'x', y: 'y'});
    view.dockManager.dock(sp, DG.DOCK_TYPE.RIGHT, null);
    await Promise.all([
      testEvent(sp.onAfterDrawScene, () => {}, () => { sp.invalidateCanvas(); }, 1000),
      awaitGrid(view.grid, 500)
    ]);

    const spBcr = sp.root.getBoundingClientRect();
    const wp = sp.worldToScreen(1, 0);
    const ev = new MouseEvent('mousemove', {
      cancelable: true, bubbles: true, view: window, button: 0,
      clientX: spBcr.left + wp.x, clientY: spBcr.top + wp.y
    });
    const spCanvas = $(sp.root).find('canvas').get()[0] as HTMLCanvasElement;
    await testEvent(fromEvent(spCanvas, 'mousemove'), () => {
      _package.logger.debug(`Test: event, currentRowIdx=${df.currentRowIdx}`);
      expect($(ui.tooltip.root).find('div table.d4-row-tooltip-table tr td canvas').length, 1);
      expect(sp.hitTest(wp.x, wp.y), 1);
    }, () => {
      spCanvas.dispatchEvent(ev);
    }, 500);
    // TODO: Any error occurred become 'Cannot read properties of null (reading 'get$columns')' because of scatter plot
    //await testEvent(sp.onAfterDrawScene, () => {}, () => { sp.invalidateCanvas(); }, 200);
    await awaitGrid(view.grid, 1000);
  }
});
