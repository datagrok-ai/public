import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {category, before, after, test, delay, awaitCheck} from '@datagrok-libraries/utils/src/test';

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

    await awaitCheck(()=> {
      const divs = Array.from(it.querySelectorAll('div'))
        .filter((it) => (it as HTMLElement).innerText === '1-methyl-4-[4-(trifluoromethyl)phenoxy]piperidine');
      return divs.length > 0;
    }, 'Incorrect name', 3000);

    ih?.click(); await delay(10);
    v.close();
    (document.querySelector('.fa-chevron-square-up') as HTMLElement)?.click();
    grok.shell.o = ui.div();
  });


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
