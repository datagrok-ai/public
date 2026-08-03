import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

import {after, before, category, delay, expect, test, testEvent} from '@datagrok-libraries/test/src/test';
import {IHelmHelper, getHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/types/monomer-library';
import {
  getUserLibSettings, setUserLibSettings,
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
    await monomerLibHelper.loadMonomerLibForTests(); // load default libraries
  });

  after(async () => {
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true); // load user settings libraries
  });

  test('inDialog', async () => {
    const helmValue = 'PEPTIDE1{[meY].A.G.T}$$$$';
    const helmInput = await ui.input.helmAsync('Macromolecule');
    helmInput.stringValue = helmValue;
    const dlg = ui.dialog('Helm Input Test')
      .add(helmInput)
      .show();
    await delay(500);

    const mol = helmInput.molValue;
    expect(mol.atoms.length, 4);
    expect(mol.bonds.length, 3);
  });

  test('renders-with-nonzero-height (root in a dialog)', async () => {
    // Regression: `ui.input.helmAsync` placed in a content-sized container
    // (no explicit height) used to collapse the editor host to 0 px, rendering
    // an empty div. Assert the editor host AND its mounted SVG are visibly tall.
    const helmInput = await ui.input.helmAsync('Macromolecule');
    helmInput.stringValue = 'PEPTIDE1{[meY].A.G.T}$$$$';
    const dlg = ui.dialog('Helm Input Height Test').add(helmInput).show();
    try {
      await delay(500);
      const hostEl = helmInput.getInput();
      expect(hostEl.clientHeight > 0, true, `editor host height must be > 0 (was ${hostEl.clientHeight})`);
      const svg = hostEl.querySelector('svg');
      expect(svg != null, true, 'an SVG must be mounted in the editor host');
      const svgH = (svg as SVGSVGElement).getBoundingClientRect().height;
      expect(svgH > 0, true, `mounted SVG height must be > 0 (was ${svgH})`);
    } finally {
      dlg.close();
    }
  });

  test('default size is 250x250', async () => {
    // `ui.input.helmAsync` with no `editorOptions` should render a compact
    // 250×250 box, NOT stretch to the form row or fall back to hwe's 400×100.
    const helmInput = await ui.input.helmAsync('Macromolecule');
    helmInput.stringValue = 'PEPTIDE1{[meY].A.G.T}$$$$';
    const dlg = ui.dialog('Helm Input Default Size').add(helmInput).show();
    try {
      await delay(500);
      const host = helmInput.getInput();
      expect(Math.abs(host.clientWidth - 250) <= 4, true, `default width ~250 (was ${host.clientWidth})`);
      expect(Math.abs(host.clientHeight - 250) <= 4, true, `default height ~250 (was ${host.clientHeight})`);
    } finally {
      dlg.close();
    }
  });

  test('editorOptions.width/height override the default size', async () => {
    const helmInput = await ui.input.helmAsync('Macromolecule', {editable: true, editorOptions: {width: 320, height: 180}});
    helmInput.stringValue = 'PEPTIDE1{[meY].A.G.T}$$$$';
    const dlg = ui.dialog('Helm Input Custom Size').add(helmInput).show();
    try {
      await delay(500);
      const host = helmInput.getInput();
      expect(Math.abs(host.clientWidth - 320) <= 4, true, `width ~320 from editorOptions (was ${host.clientWidth})`);
      expect(Math.abs(host.clientHeight - 180) <= 4, true, `height ~180 from editorOptions (was ${host.clientHeight})`);
    } finally {
      dlg.close();
    }
  });

  test('editor canvas sizes from getInput() box (consumer sets height on the host)', async () => {
    // The real consumer pattern (e.g. libraries/db-explorer/renderer.ts) sets an
    // explicit size directly on `getInput()` (the editor host). The editor
    // canvas must follow the HOST's own box — the old code derived the size from
    // `this.root` and painted a 0-height SVG when a consumer sized only the host
    // (the reported "height 0 / width 105" bug). Mount via the input so layout
    // is real, then size the host and assert the SVG tracks it.
    const helmInput = await ui.input.helmAsync('Macromolecule');
    helmInput.stringValue = 'PEPTIDE1{[meY].A.G.T}$$$$';
    const dlg = ui.dialog('Helm getInput() Size Test').add(helmInput).show();
    try {
      await delay(300);
      const host = helmInput.getInput();
      // Consumer sizes the host element directly (the db-explorer pattern).
      host.style.setProperty('width', '360px', 'important');
      host.style.setProperty('height', '260px', 'important');
      await delay(500); // ResizeObserver → editor.setSize → re-render

      expect(host.clientWidth, 360, `host width (was ${host.clientWidth})`);
      expect(host.clientHeight, 260, `host height (was ${host.clientHeight})`);
      const svg = host.querySelector('svg') as SVGSVGElement | null;
      expect(svg != null, true, 'an SVG must be mounted');
      const r = svg!.getBoundingClientRect();
      // The canvas tracks the host box (NOT root): both dimensions filled.
      expect(r.height >= host.clientHeight - 8, true,
        `SVG height should track the host's ${host.clientHeight} (was ${r.height})`);
      expect(r.width >= host.clientWidth - 8, true,
        `SVG width should track the host's ${host.clientWidth} (was ${r.width})`);
    } finally {
      dlg.close();
    }
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
  const tmBcr = grok.shell.v ? grok.shell.v.root.getBoundingClientRect() : null;
  const dlg = ui.dialog('Helm Input Test')
    .add(helmInput)
    .show({...(tmBcr ? {x: tmBcr.right + 20, y: tmBcr.top + 20} : {})});
  try {
    const dialogInputEl = $(dlg.root)
      .find('div.d4-dialog-contents > div.ui-input-helm > div.ui-input-editor').get(0);
    expect(helmInput.getInput(), dialogInputEl);

    const mon = {i: 0, elem: 'meY'};
    const a = helmInput.molValue.atoms[mon.i];
    const iEl = helmInput.getInput();
    const iBcr = iEl.getBoundingClientRect();
    const ev = new MouseEvent('mousemove', {
      cancelable: true, bubbles: true, view: window, button: 0,
      clientX: iBcr.left + a.p.x, clientY: iBcr.top + a.p.y
    });
    await testEvent(ui.tooltip.onTooltipShown, () => { /* tooltip shown*/ },
      () => {
        iEl.dispatchEvent(ev);
        // do not await after event dispatch before test tooltip content
        // await leads to tooltip changed (with TestManager tree node tooltip)
      }, 100, 'timeout await ui.tooltip.onTooltipShown');
    const tooltipEl = $(document).find('body > div.d4-tooltip').get(0) as HTMLElement;
    const elemDiv = $(tooltipEl).find('div > div.d4-flex-row.ui-div > div:nth-child(1)').get(0);
    expect(tooltipEl.style.display != 'none', true);
    expect(elemDiv?.innerText, mon.elem);
  } finally {
    await delay(50);
    ui.tooltip.hide();
  }
}
