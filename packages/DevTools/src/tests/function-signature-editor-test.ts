import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import { category, before, test, expect, waitUntil, after } from '@datagrok-libraries/test/src/test';

category('FSE', () => {
  const script: DG.Script = DG.Script.create(`#name: Template
#description: Calculates number of cells in the table
#language: python
#sample: cars.csv
#input: dataframe table [Data table]
#output: int count [Number of cells in table]

count = table.shape[0] * table.shape[1]`);

  before(async () => {
    grok.shell.windows.simpleMode = false;
  });

  test('exist', async () => { 
    let b: boolean;
    for (let i = 0; i < 5; i++) {
      b = true;
      try {
        const v = DG.ScriptView.create(script);
        grok.shell.addView(v);
        await waitUntil(() => grok.shell.v != null, '', 1000);
        await waitUntil(() => {
          return !!(grok.shell.v.getRibbonPanels().flat().find(
            (elem) => !!elem.querySelector('.fa-magic'),
          )?.lastChild as HTMLElement);
        }, 'Failed: FSE button does not exist', 1500);
        break;
      } catch (e) {
        b = false;
      }
    }
    if (!b) throw new Error('Failed: FSE button does not exist');
  }, {owner:'oserhiienko@datagrok.ai'});

  test('open', async () => {
    const v = DG.ScriptView.create(script);
    grok.shell.addView(v);
    await waitUntil(() => grok.shell.v != null, '', 2000);
    let fseButton: HTMLElement;
    await waitUntil(() => {
      fseButton = grok.shell.v.getRibbonPanels().flat().find(
        (elem) => !!elem.querySelector('.fa-magic'),
      )?.lastChild as HTMLElement;
      return !!fseButton;
    }, 'Failed: FSE button does not exist', 2000);
    fseButton.click();
    await waitUntil(() => !!grok.shell.v.root.querySelector('[name="PARAMETERS"]'),
      'FSE button did not open view', 2000);
    const currentView = grok.shell.v;
    const allTabsExist =
      !!currentView.root.querySelector('[name="PROPERTIES"]') &&
      !!currentView.root.querySelector('[name="PARAMETERS"]') &&
      !!currentView.root.querySelector('[name="CODE"]') &&
      !!currentView.root.querySelector('[name="UI"]');
    expect(allTabsExist, true);
  }, {owner:'oserhiienko@datagrok.ai'});

  test('close', async () => {
    const v = DG.ScriptView.create(script);
    grok.shell.addView(v);
    await waitUntil(() => grok.shell.v != null, '', 2000);
    let fseButton: HTMLElement;
    await waitUntil(() => {
      fseButton = grok.shell.v.getRibbonPanels().flat().find(
        (elem) => !!elem.querySelector('.fa-magic'),
      )?.lastChild as HTMLElement;
      return !!fseButton;
    }, 'Failed: FSE button does not exist', 2000);
    fseButton.click();
    await waitUntil(() => !!grok.shell.v.root.querySelector('[name="PARAMETERS"]'),
      'FSE button did not open view', 2000);
    const currentView = grok.shell.v;
    let editorButton: HTMLElement;
    await waitUntil(() => {
      editorButton = grok.shell.v.getRibbonPanels().flat().find(
        (elem) => !!elem.querySelector('.fa-code'),
      )?.lastChild as HTMLElement;
      return !!editorButton;
    }, 'Failed: FSE button does not exist', 2000);
    editorButton.click();
    await waitUntil(() => !!currentView.root.querySelector('.CodeMirror.cm-s-default'),
      'Code button did not open view', 2000);
  }, {owner:'oserhiienko@datagrok.ai'});
});
