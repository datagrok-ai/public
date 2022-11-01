import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as bio from '@datagrok-libraries/bio';

import {after, before, category, test, expect, expectArray, expectObject} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import {newickToDf} from '../utils';
import {PhylocanvasGlViewer} from '../viewers/phylocanvas-gl-viewer';

category('phylocanvasGlViewer', () => {

  let viewList: DG.ViewBase[];
  let dfList: DG.DataFrame[];
  let currentView: DG.ViewBase;

  before(async () => {
    viewList = [];
    dfList = [];
    currentView = grok.shell.v;
  });

  after(async () => {
    viewList.forEach((v) => { v.close(); });
    dfList.forEach((df) => { grok.shell.closeTable(df); });
    grok.shell.v = currentView;

  });

  test('open', async () => {
    await _testOpen();
  });

  async function _testOpen() {
    const newickStr = await _package.files.readAsText(`data/tree95.nwk`);
    const treeDf: DG.DataFrame = newickToDf(newickStr, 'tree95');

    const tv = grok.shell.addTableView(treeDf, DG.DOCK_TYPE.FILL);

    const viewer: bio.IPhylocanvasGlViewer = (await treeDf.plot.fromType('PhylocanvasGl', {
      interactive: true,
      alignLabels: true,
      showLabels: true,
      showLeafLabels: true,
      padding: 0,
      treeToCanvasRatio: 1,
    })) as unknown as bio.IPhylocanvasGlViewer;

    const viewerOnAfterRenderPromise = new Promise<void>((resolve, reject) => {
      try {
        viewer.onAfterRender.subscribe(({gl}) => {
          resolve();
        });
      } catch (err) {
        reject(err);
      }
    });
    const treeDn = tv.dockManager.dock(viewer as unknown as DG.Viewer, DG.DOCK_TYPE.LEFT);
    await viewerOnAfterRenderPromise;
    let k = 11;
  }
});