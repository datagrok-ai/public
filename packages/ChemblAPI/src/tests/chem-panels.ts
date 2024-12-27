import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {category, before, after, expect, test, delay, awaitCheck} from '@datagrok-libraries/utils/src/test';

category('UI info panel', () => {
  let v: DG.TableView;
  let smiles: DG.DataFrame;

  before(async () => {
    grok.shell.closeAll();
    grok.shell.windows.showProperties = true;
  });


  test('identifiers', async () => {
    smiles = grok.data.demo.molecules(20);
    v = grok.shell.addTableView(smiles);
    await grok.data.detectSemanticTypes(smiles);
    await awaitCheck(() => document.querySelector('canvas') !== null, 'cannot load table', 3000);
    const pp = document.querySelector('.grok-prop-panel') as HTMLElement;
    await awaitPanel(pp, 'Structure');
    (document.querySelector('.fa-chevron-square-down') as HTMLElement)?.click();
    await awaitPanel(pp, 'Identifiers', 2000);
    const ih = Array.from(pp.querySelectorAll('div.d4-accordion-pane-header'))
      .find((el) => el.textContent === 'Identifiers') as HTMLElement;
    await delay(200);
    if (!ih?.classList.contains('expanded')) ih?.click();
    await awaitCheck(() => (ih.nextSibling as HTMLElement)
      .querySelector('table') !== null, 'cannot load Identifiers', 15000);
    const it = ih.nextSibling as HTMLElement;

    const links = it.querySelectorAll('.ui-link');
    for (const i of ['SCHEMBL5536145', '18722989', 'CHEMBL2262190']) {
      let found = false;
      innerLoop:
      for (let j = 0; j < links.length; j++) {
        if (links[j].textContent === i) {
          found = true;
          break innerLoop;
        }
      }
      expect(found, true, `${i} identifier not found`);
    }

    ih?.click(); await delay(10);
    v.close();
    (document.querySelector('.fa-chevron-square-up') as HTMLElement)?.click();
    grok.shell.o = ui.div();
  });

  test('map identifiers', async () => {
    smiles = grok.data.demo.molecules(20);
    v = grok.shell.addTableView(smiles);
    await grok.data.detectSemanticTypes(smiles);
    await awaitCheck(() => document.querySelector('canvas') !== null, 'cannot load table', 3000);

    await callDialog();
    setDialogInputValue('Chem Map Identifiers', 'To Source', 'mcule');
    let okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
    okButton!.click();
    await awaitCheck(() => grok.shell.t.columns.contains('mcule'), 'cannot find mcule column', 15000);

    await callDialog();
    setDialogInputValue('Chem Map Identifiers', 'To Source', 'chembl');
    okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
    okButton!.click();
    await awaitCheck(() => grok.shell.t.columns.contains('chembl'), 'cannot find chembl column', 15000);

    await callDialog();
    setDialogInputValue('Chem Map Identifiers', 'To Source', 'pubchem');
    okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
    okButton!.click();
    await awaitCheck(() => grok.shell.t.columns.contains('pubchem'), 'cannot find pubchem column', 15000);

    v.close();
    grok.shell.o = ui.div();

    async function callDialog() {
      grok.shell.topMenu.find('Chem').group('Calculate').find('Map Identifiers...').click();
      await awaitCheck(() => DG.Dialog.getOpenDialogs().length > 0, 'cannot find Chem Map Identifiers dialog', 1000);
    }
  }, {timeout: 60000});


  after(async () => {
    grok.shell.closeAll();
  });
});

async function awaitPanel(pp: HTMLElement, name: string, ms: number = 5000): Promise<void> {
  await awaitCheck(() => {
    return Array.from(pp?.querySelectorAll('div.d4-accordion-pane-header'))
      .find((el) => el.textContent === name) !== undefined;
  }, `cannot find ${name} property`, ms);
}

export function returnDialog(dialogTitle:string):DG.Dialog | undefined {
  let dialog: DG.Dialog | undefined;
  for (let i = 0; i < DG.Dialog.getOpenDialogs().length; i++) {
    if (DG.Dialog.getOpenDialogs()[i].title == dialogTitle) {
      dialog = DG.Dialog.getOpenDialogs()[i];
      return dialog;
    }
  }
}

export function setDialogInputValue(dialogTitle:string, inputName:string,
  value:string | number | boolean | DG.Column | DG.Column[]):void {
      returnDialog(dialogTitle)!.input(inputName).value = value;
}
