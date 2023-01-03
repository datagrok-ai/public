import {after, before, category, test, awaitCheck} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {checkHTMLElement} from './utils';

category('UI: Tables', () => {
  let v: DG.View;
  const df = grok.data.demo.demog(100);

  const arr = [
    {key: 'First', value: true},
    {key: 'Second', value: false},
    {key: 'Third', value: false},
  ];

  const table = ui.table(arr, (item, idx) => [item.key, item.value]);
  const tableFromMap = ui.tableFromMap(arr);
  const HTMLTable = DG.HtmlTable.create(arr, (item: any, idx: any) => [item.key, item.value]);

  before(async () => {
    v = grok.shell.newView('');
    grok.shell.windows.showProperties = true;
  });

  test('table.root', async () => {
    checkHTMLElement('Table', table, v, 'table');
  });

  test('tableFromMap.root', async () => {
    checkHTMLElement('tableFromMap', tableFromMap, v, 'table');
  });

  test('HTMLTable.root', async () => {
    checkHTMLElement('HTMLTable', HTMLTable.root, v, 'table');
  });

  test('normalize (min-max)', async () => {
    await normalize();
  });

  test('normalize (z-scores)', async () => {
    await normalize(true);
  });

  after(async () => {
    v.close();
    grok.shell.closeAll();
  });

  async function normalize(changeMethod: boolean = false) {
    const ageDF = DG.DataFrame.fromColumns([df.getCol('age')]);
    const tv = grok.shell.addTableView(ageDF);
    try {
      const pp = document.querySelector('.grok-prop-panel') as HTMLElement;
      await awaitCheck(() => {
        if (pp.querySelector('.ui-label') !== null)
          return (pp.querySelector('.ui-label') as HTMLElement).innerText === 'age';
        return false;
      }, 'cannot load property panel of age column', 5000);
      const actions = Array.from(pp.querySelectorAll('div.d4-accordion-pane-header'))
        .find((el) => el.textContent === 'Actions') as HTMLElement;
      if (!actions.classList.contains('expanded')) await actions.click();
      const normalize = Array.from(pp.querySelectorAll('.d4-link-action'))
        .find((el) => el.textContent === 'Normalize...') as HTMLElement;
      normalize.click();
      await awaitCheck(() => DG.Dialog.getOpenDialogs().length > 0);
      const dialog = DG.Dialog.getOpenDialogs()[0];
      if (changeMethod) {
        const method = dialog.inputs[2];
        method.value = method.property.choices[1];
      }
      (dialog.root.querySelector('.ui-btn.ui-btn-ok.enabled') as HTMLElement).click();
      await awaitCheck(() => ageDF.getCol('age').get(3) === (changeMethod ? 0.6055613160133362 : 0.7254902124404907),
        'normalize does not work properly', 1000);
    } finally {
      tv.close();
      grok.shell.closeTable(ageDF);
    }
  }
});
