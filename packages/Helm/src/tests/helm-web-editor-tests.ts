import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, awaitCheck, before, category, delay, expect, test, expectArray} from '@datagrok-libraries/test/src/test';
import {IHelmHelper, getHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/types/monomer-library';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {
  getUserLibSettings, setUserLibSettings
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';

import {initHelmMainPackage} from './utils';

category('helm-web-editor', async () => {
  let monomerLibHelper: IMonomerLibHelper;
  let userLibSettings: UserLibSettings;

  let hh: IHelmHelper;

  before(async () => {
    await initHelmMainPackage();
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();

    await monomerLibHelper.loadMonomerLibForTests();

    hh = await getHelmHelper();
  });

  after(async () => {
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true);
  });

  test('helm-web-editor', async () => {
    const editor = hh.createHelmWebEditor();
  });

  test('web-editor-app', async () => {
    const helmStr: string = 'PEPTIDE1{I.H.A.N.T.Thr_PO3H2.Aca.D-Tyr_Et}$$$$';

    const appHost = ui.div([],
      {style: {width: '700px', height: '432px'}});
    const app = hh.createWebEditorApp(appHost, helmStr);

    const dlg = ui.dialog().add(appHost).show();
    try {

    } finally {
      dlg.close();
    }
  });

  test('a dialog closes after its OK handler', async () => {
    // The editor's dialogs read the editor on OK and destroy it on close: OK has to come first
    let okRan = false;
    let okBeforeClose: boolean | null = null;
    const dlg = ui.dialog().add(ui.div()).onOK(() => { okRan = true; });
    const closeSub = dlg.onClose.subscribe(() => {
      closeSub.unsubscribe();
      okBeforeClose = okRan;
    });
    dlg.show();
    (dlg.root.querySelector('.ui-btn-ok') as HTMLElement).click();
    await awaitCheck(() => okBeforeClose !== null, 'the dialog did not close', 3000);
    expect(okBeforeClose, true);
  });

  test('cell editor: Escape takes the monomer that follows the mouse away with the dialog', async () => {
    // A palette tile clicked, then Escape: the dialog closed and the monomer that follows the mouse stayed on the page
    // (4 November). The dialog never destroyed the editor, and its canvas left the page without a mouseleave.
    const df = DG.DataFrame.fromCsv('helm\nPEPTIDE1{A.C.D}$$$$\n');
    await grok.data.detectSemanticTypes(df);
    const tv = grok.shell.addTableView(df);
    const appRootOnPage = (): HTMLElement | null =>
      document.querySelector('.d4-dialog [data-testid="app-root"]') as HTMLElement | null;
    try {
      // the cell editor a double click opens
      await grok.functions.call('Helm:editMoleculeCell', {cell: tv.grid.cell('helm', 0)});
      await awaitCheck(() => appRootOnPage() !== null, 'the editor did not open', 5000);
      const appRoot = appRootOnPage()!;
      (appRoot.querySelector('[data-testid^="palette-tile-"]') as HTMLElement).click();
      const canvasHost = appRoot.querySelector('[data-testid="canvas-host"]') as HTMLElement;
      const r = canvasHost.getBoundingClientRect();
      // a corner of the canvas, off the three monomers
      canvasHost.dispatchEvent(
        new MouseEvent('mousemove', {clientX: r.right - 10, clientY: r.bottom - 10, bubbles: true}));
      const preview = document.querySelector('body > [data-testid="next-add-cursor"]') as HTMLElement | null;
      expect(preview?.style.display, 'block');

      const esc = new KeyboardEvent('keydown', {key: 'Escape', code: 'Escape', bubbles: true, cancelable: true});
      Object.defineProperty(esc, 'keyCode', {get: () => 27});
      Object.defineProperty(esc, 'which', {get: () => 27});
      (document.activeElement ?? document.body).dispatchEvent(esc);
      await awaitCheck(() => appRootOnPage() === null, 'Escape did not close the dialog', 3000);
      expect(preview!.style.display, 'none');
    } finally {
      for (const d of DG.Dialog.getOpenDialogs()) d.close();
      tv.close();
    }
  });
});
