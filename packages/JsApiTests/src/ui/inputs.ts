import { after, before, category, delay, expect, test } from "../test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { caption, enabled, HTMLElement } from './utils';


category('UI: Inputs', () => {
  let v: DG.View;
  let t = grok.data.testData('demog', 100);
  let tables = grok.shell.tables;
  const inputs = {
    'stringInput': ui.stringInput('', ''),
    'intInput': ui.intInput('', 0),
    'floatInput': ui.floatInput('', 0.00),
    'boolInput': ui.boolInput('', true),
    'switchInput': ui.switchInput('', true),
    'choiceInput': ui.choiceInput('', '', []),
    'multiChoiceInput': ui.multiChoiceInput('', [], []),
    'dateInput': ui.dateInput('', DG.DateTime.fromDate(new Date(2000, 1, 1))),
    'textInput': ui.textInput('', ''),
    'searchInput': ui.searchInput('', ''),
    'columnInput': ui.columnInput('', t, t.col('age')),
    'tableInput': ui.tableInput('', tables[0], tables, (table: any) => grok.shell.info(table.name))
  };

  before(async () => {
    v = grok.shell.newView('');
  });

  test('input.root', async () => {
    for (const [key, value] of Object.entries(inputs)) {
      HTMLElement(key, value, v, '.ui-input-root');
    }
  })

  test('input.input', async () => {
    for (const [key, value] of Object.entries(inputs)) {
      HTMLElement(key, value, v, '.ui-input-editor');
    }
  })

  test('input.captionLabel', async () => {
    for (const [key, value] of Object.entries(inputs)) {
      HTMLElement(key, value, v, '.ui-input-label');
    }
  })

  test('input.caption', async () => {
    for (const [key, value] of Object.entries(inputs)) {
      caption(key, value, v, '.ui-input-label');
    }
  })

  test('input.stringValue', async () => {
    for (const [key, value] of Object.entries(inputs)) {
      stringValue(key, value, '.ui-input-editor');
    }
  })

  test('input.ebabled', async () => {
    for (const [key, value] of Object.entries(inputs)) {
      enabled(key, value, v, '.ui-input-root');
    }
  })

  after(async () => {
    v.close();
    grok.shell.closeAll();
  });

  function stringValue(name: string, input: DG.InputBase, selector: string): void {
    v.append(input.root);
    let value: string;
    try {
      value = (<HTMLInputElement>v.root.querySelector(selector)).value;
      expect(input.stringValue, value)
    }
    catch (x) {
      throw name + ': ' + x;
    }
    finally {
      input.root.remove();
    }
  }

});