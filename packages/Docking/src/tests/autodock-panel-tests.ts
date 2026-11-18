import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, expect, test, waitUntil} from '@datagrok-libraries/test/src/test';

import {_package} from '../package-test';

category('AutoDock panel', () => {
  test('receptor with the pose', async () => {
    const df = DG.DataFrame.fromCsv(await _package.files.readAsText('demo_files/docking_demo.csv'));
    const poses = df.getCol('AutoDock poses');
    poses.semType = DG.SEMTYPE.MOLECULE3D;
    poses.meta.units = 'pdb';

    const widget: DG.Widget = await grok.functions.call('Docking:getAutodockSingle', {
      molecule: DG.SemanticValue.fromTableCell(df.cell(0, poses.name)), showProperties: false, table: df});
    const dialog = ui.dialog('AutoDock').add(widget.root).show();
    try {
      const viewer = DG.Widget.find(widget.root.querySelector('.bsv-container-info-panel')!) as any;
      await waitUntil(() => !viewer.isRenderPending, 'the Biostructure viewer did not finish loading', 30000);
      expect(viewer.dataEff?.options?.name, 'BACE1', 'the viewer shows the pose instead of the receptor');
      expect(viewer.ligandColumnName, poses.name);
    } finally {
      dialog.close();
    }
  }, {timeout: 60000});
}, {owner: 'oserhiienko@datagrok.ai'});
