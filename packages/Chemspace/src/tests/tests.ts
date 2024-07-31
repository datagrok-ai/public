import * as grok from 'datagrok-api/grok';
import {test, category, awaitCheck} from '@datagrok-libraries/utils/src/test';
import {app, pricesPanel, samplesPanel} from '../package';

category('Chemspace', () => {
  const mol = 'Oc1ccccc1';

  test('Prices panel', async () => {
    const widget = await pricesPanel('CSCS00102564062');
    await awaitCheck(() => widget.root.getElementsByClassName('chemspace-prices-grid').length > 0,
      'prices panel hasn\'t been created', 30000);
  }, {skipReason: 'Requires API key'});

  test('Samples panel', async () => {
    const widget = await samplesPanel(mol);
    const similarPaneHeader = widget.root.querySelector('[name="div-section--Similar"]') as HTMLElement;
    similarPaneHeader.click();
    const similarPane = widget.root.querySelector('[d4-title="Similar"]') as HTMLElement;
    await awaitCheck(() => similarPane.getElementsByClassName('chem-mol-box').length === 9,
      `Expected 9 similar molecules, got ${similarPane.getElementsByClassName('chem-mol-box').length}`, 30000);
    const subPaneHeader = widget.root.querySelector('[name="div-section--Substructure"]') as HTMLElement;
    subPaneHeader.click();
    const subPane = widget.root.querySelector('[d4-title="Substructure"]') as HTMLElement;
    await awaitCheck(() => subPane.getElementsByClassName('chem-mol-box').length === 10,
      `Expected 9 molecules with substructure, got ${subPane.getElementsByClassName('chem-mol-box').length}`, 30000);

  }, { skipReason: 'Requires API key' , timeout: 60000});

  test('App', async () => {
    app();
    await awaitCheck(() => grok.shell.tv.dataFrame.rowCount === 10,
    `Search hasn't been completed`, 30000);
  }, {skipReason: 'Requires API key'});
});
