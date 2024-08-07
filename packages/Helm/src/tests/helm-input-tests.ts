import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

import {after, before, category, delay, expect, test, timeout} from '@datagrok-libraries/utils/src/test';
import {IHelmHelper, getHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {
  getUserLibSettings, setUserLibSettings, setUserLibSettingsForTests
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';

import type {HelmInput} from '../widgets/helm-input';

category('HelmInput', () => {
  let helmHelper: IHelmHelper;

  let monomerLibHelper: IMonomerLibHelper;
  /** Backup actual user's monomer libraries settings */
  let userLibSettings: UserLibSettings;

  before(async () => {
    helmHelper = await getHelmHelper();

    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();

    // Test 'helm' requires default monomer library loaded
    await setUserLibSettingsForTests();
    await monomerLibHelper.awaitLoaded();
    await monomerLibHelper.loadLibraries(true); // load default libraries
  });

  after(async () => {
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadLibraries(true); // load user settings libraries
  });

  test('inDialog', async () => {
    const helmValue = 'PEPTIDE1{[meY].A.G.T}$$$$';
    const helmInput = await ui.input.helmAsync('Macromolecule');
    helmInput.stringValue = helmValue;
    const dlg = ui.dialog('Helm Input Test')
      .add(helmInput)
      .show();
    await delay(500);

    const mol = helmInput.value;
    expect(mol.atoms.length, 4);
    expect(mol.bonds.length, 3);
  });

  test('tooltip', async () => {
    await _testTooltipOnHelmInput();
  }, {timeout: 300000});
});

async function _testTooltipOnHelmInput(): Promise<void> {
  const helmValue = 'PEPTIDE1{[meY].A.G.T}$$$$';
  const helmInput = (await ui.input.helmAsync('Macromolecule')) as HelmInput;
  helmInput.stringValue = helmValue;
  // Show dialog to the right of the TestManager test tree
  // to prevent interfering mousemove position for tooltip.
  const tmBcr = grok.shell.v.root.getBoundingClientRect();
  const dlg = ui.dialog('Helm Input Test')
    .add(helmInput)
    .show({x: tmBcr.right + 20, y: tmBcr.top + 20});
  try {
    const dialogInputEl = $(dlg.root)
      .find('div.d4-dialog-contents > div.ui-input-helm > div.ui-input-editor').get(0);
    expect(helmInput.getInput(), dialogInputEl);


    const mon = {i: 0, elem: 'meY'};
    const a = helmInput.value.atoms[mon.i];
    const iEl = helmInput.getInput();
    const iBcr = iEl.getBoundingClientRect();
    const ev = new MouseEvent('mousemove', {
      cancelable: true, bubbles: true, view: window, button: 0,
      clientX: iBcr.left + a.p.x, clientY: iBcr.top + a.p.y
    });
    iEl.dispatchEvent(ev);
    // do not await anything after event dispatch before test tooltip content
    // await leads to tooltip changed (with TestManager tree node tooltip)

    const tooltipEl = $(document).find('body > div.d4-tooltip').get(0) as HTMLElement;
    const elemDiv = $(tooltipEl).find('div > div.d4-flex-row.ui-div > div:nth-child(1)').get(0);
    expect(tooltipEl.style.display != 'none', true);
    expect(elemDiv?.innerText, mon.elem);
  } finally {
    await delay(50);
    ui.tooltip.hide();
  }
}
