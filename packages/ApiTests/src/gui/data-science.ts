import {category, test, awaitCheck} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


category('Data Science', () => {
  let tv: DG.TableView;

  test('SPE 2D', async () => {
    await spe('2D');
  });

  test('SPE 3D', async () => {
    await spe('3D');
  });

  async function spe(dimensions: string) {
    const smiles = grok.data.demo.molecules(100);
    tv = grok.shell.addTableView(smiles);
    try {
      await awaitCheck(() => document.querySelector('canvas') !== null, 'cannot load table', 3000);
      grok.shell.topMenu.find('Tools').find('Data Science').find('Stochastic Proximity Embedding...').click();
      await awaitCheck(() => DG.Dialog.getOpenDialogs().length === 1,
        'cannot find Stochastic Proximity Embedding dialog', 1000);
      const dialog = DG.Dialog.getOpenDialogs()[0];
      dialog.input('Columns').value = [smiles.col('smiles')];
      dialog.input('Dimensions').value = dimensions;
      dialog.input('Columns').input.dispatchEvent(new MouseEvent('click'));
      await awaitCheck(() => DG.Dialog.getOpenDialogs().length === 2, 'cannot find Select columns dialog', 1000);
      (DG.Dialog.getOpenDialogs()[1].root.querySelector('.ui-btn.ui-btn-ok.enabled') as HTMLElement).click();
      (dialog.root.querySelector('.ui-btn.ui-btn-ok.enabled') as HTMLElement).click();
      await awaitCheck(() => smiles.col('x') !== null, 'cannot find x column', 3000);
      await awaitCheck(() => smiles.col('y') !== null, 'cannot find y column', 1000);
      if (dimensions === '3D') await awaitCheck(() => smiles.col('z') !== null, 'cannot find z column', 1000);
    } finally {
      tv.close();
      grok.shell.closeTable(smiles);
      DG.Dialog.getOpenDialogs().forEach((d) => d.close());
    }
  }
});
