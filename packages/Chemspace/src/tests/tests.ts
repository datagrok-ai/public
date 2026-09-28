import * as grok from 'datagrok-api/grok';
import {test, category, awaitCheck, expect} from '@datagrok-libraries/test/src/test';

category('Chemspace', () => {
  const mol = 'Oc1ccccc1';

  async function expectMolCount(pane: HTMLElement, expected: number): Promise<void> {
    const count = (): number => pane.getElementsByClassName('chem-mol-box').length;
    await awaitCheck(() => count() >= expected, `Molecules haven't been rendered`, 30000);
    expect(count(), expected);
  }

  test('Prices panel', async () => {
    const widget = await grok.functions.call('Chemspace:pricesPanel', {'id': 'CSSS00102643788'});
    await awaitCheck(() => widget.root.getElementsByClassName('chemspace-prices-grid').length > 0,
      'prices panel hasn\'t been created', 30000);
  });

  test('Samples panel', async () => {
    const widget = await grok.functions.call('Chemspace:samplesPanel', {'smiles': mol});
    const similarPaneHeader = widget.root.querySelector('[name="div-section--Similar"]') as HTMLElement;
    similarPaneHeader.click();
    // Chemspace returns 10 hits; the 14C-labelled phenol scores below the panel's similarity cutoff
    await expectMolCount(widget.root.querySelector('[d4-title="Similar"]') as HTMLElement, 9);
    const subPaneHeader = widget.root.querySelector('[name="div-section--Substructure"]') as HTMLElement;
    subPaneHeader.click();
    await expectMolCount(widget.root.querySelector('[d4-title="Substructure"]') as HTMLElement, 10);
  }, {timeout: 60000});

  test('App', async () => {
    await grok.functions.call('Chemspace:app');
    await awaitCheck(() => grok.shell.tv.dataFrame.rowCount > 0,
      `Search hasn't been completed`, 30000);
  }, {timeout: 60000});
});
