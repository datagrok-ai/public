import { after, before, category, expect, expectArray, test } from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { caption, enabled, checkHTMLElement, stringValue } from './utils';
import dayjs from 'dayjs';


category('UI: Inputs', () => {
  let v: DG.View;
  const t = grok.data.testData('demog', 100);
  const tables = grok.shell.tables;
  let inputs: { [key: string]: DG.InputBase };

  before(async () => {
    inputs = {
      'stringInput': ui.input.string('', { value: '' }),
      'intInput': ui.input.int('', { value: 0 }),
      'floatInput': ui.input.float('', { value: 0.00 }),
      'boolInput': ui.input.bool('', { value: true }),
      'switchInput': ui.input.toggle('', { value: true }),
      'choiceInput': ui.input.choice('', { value: '1', items: ['1', '2', '3'] }),
      'multiChoiceInput': ui.input.multiChoice('', { value: [], items: [] }),
      // 'dateInput': ui.input.date('', {value: dayjs('2001-01-01')}),
      'textInput': ui.input.textArea('', { value: '' }),
      'searchInput': ui.input.search('', { value: '' }),
      'columnInput': ui.input.column('', { table: t, value: t.col('age')! }),
      'columnsInput': ui.input.columns('', { table: t, onValueChanged: () => null }),
      'tableInput': ui.input.table('', { value: tables[0], items: tables }),
      'colorInput': ui.input.color('', { value: '#ff0000' }),
    };
    v = grok.shell.newView('');

    for (let inputName of Object.keys(inputs)) {
      v.append(inputs[inputName]);
    }
    grok.shell.windows.showContextPanel = true;
  });

  test('input.root', async () => {
    for (const [key, value] of Object.entries(inputs))
      checkHTMLElement(key, value.root, v, '.ui-input-root');
  });

  test('input.input', async () => {
    for (const [key, value] of Object.entries(inputs))
      checkHTMLElement(key, value.root, v, '.ui-input-editor');
  });

  test('input.captionLabel', async () => {
    for (const [key, value] of Object.entries(inputs))
      checkHTMLElement(key, value.root, v, '.ui-input-label');
  });

  test('input.caption', async () => {
    for (const [key, value] of Object.entries(inputs))
      caption(key, value, v, '.ui-input-label');
  });

  test('input.i', async () => {
    for (const [key, value] of Object.entries(inputs))
      stringValue(key, value, '.ui-input-editor', v);
  });

  test('input.enabled', async () => {
    for (const [key, value] of Object.entries(inputs))
      enabled(key, value, v, '.ui-input-root');
  });

  test('input.readOnly', async () => {
    for (const [key, value] of Object.entries(inputs))
      readonly(key, value, 'input[readonly], input[disabled]');
  });

  test('input.onChanged', async () => {
    for (const [key, value] of Object.entries(inputs)) {
      switch (key) {
        case 'stringInput':
        case 'textInput':
        case 'searchInput':
          onChanged(key, value, 'test');
          break;
        case 'intInput':
        case 'floatInput':
          onChanged(key, value, 100);
          break;
        case 'boolInput':
          onChanged(key, value, false);
          break;
        case 'dateInput':
          onChanged(key, value, dayjs('2022-01-01'));
          break;
        case 'choiceInput':
          onChanged(key, value, '2');
          break;
        case 'columnInput':
          onChanged(key, value, t.col('height'));
          break;
        case 'columnsInput':
          onChanged(key, value, [t.col('height'), t.col('weight')]);
          break;
        case 'tableInput':
          grok.shell.addTableView(grok.data.demo.demog());
          grok.shell.addTableView(grok.data.demo.randomWalk());
          onChanged(key, value, [t.col('height'), t.col('weight')]);
          break;
      }
    }
  });

  test('floatInput', async () => {
    const t = ui.input.float('Label', { value: 0.003567 });
    t.format = '0.0000';
    const v = grok.shell.newView('Test', [t]);
    const input = t.input as HTMLInputElement;
    expect(input.value, '0.0036');
    t.input.dispatchEvent(new Event('focus'));
    expect(input.value, '0.0036');
    t.value = 0.3;
    expect(input.value, '0.3000');
    v.close();
  });

  after(async () => {
    grok.shell.closeAll();
  });

  function readonly(name: string, input: DG.InputBase, selector: string): void {
    v.append(input.root);
    let value: boolean = false;
    input.readOnly = true;
    try {
      if (input.input.tagName == 'INPUT') {
        if (v.root.querySelector(selector) || name === 'switchInput')
          value = true;
        expect(input.readOnly, value);
      }
    } catch (x) {
      throw new Error(name + ': ' + x);
    } finally {
      input.readOnly = false;
      input.root.remove();
    }
  }

  function onChanged(name: string, input: DG.InputBase, newVal: any): void {
    v.append(input.root);
    const value = input.value;

    let changed = false;
    input.value = newVal;
    input.onChanged.subscribe(() => {
      changed = true;
    });

    input.fireChanged();

    if (!changed)
      throw new Error(`"${name}": OnChange value error`);

    input.value = value;
    input.fireChanged();
    input.root.remove();
  }
}, {clear: false, owner: 'dkovalyov@datagrok.ai'});

category('UI: Choice input', () => {
  test('nullable', async () => {
    const t = ui.input.choice('test', { value: '1', items: ['1', '2'], nullable: true });

    const view = grok.shell.newView();
    view.append(t);

    const selector = view.root.querySelector('select.ui-input-editor') as HTMLSelectElement;
    expect(selector.item(0)?.textContent, '');
    expect(selector.item(1)?.textContent, '1');
    expect(selector.item(2)?.textContent, '2');

    expectArray(t.items, [null, '1', '2']);
  });

  test('non-nullable', async () => {
    const t = ui.input.choice('test', { value: '1', items: ['1', '2'], nullable: false });

    const view = grok.shell.newView();
    view.append(t);

    const selector = view.root.querySelector('select.ui-input-editor') as HTMLSelectElement;
    expect(selector.item(0)?.textContent, '1');
    expect(selector.item(1)?.textContent, '2');

    expectArray(t.items, ['1', '2']);
  });

  test('fromFunction', async () => {
    const view = grok.shell.newView();
    const input = ui.input.choice('Sex', { value: 'Male', items: ['Male', 'Female'] });
    view.root.appendChild(input.root);

    input.value = null;
    expect(input.value == null);
    input.value = 'Male';
    expect(input.value, 'Male');
    input.value = 'Female';
    expect(input.value, 'Female');
    // input.value = '';
    // expect(input.value, '');

    const select = document.querySelector('select.ui-input-editor') as HTMLSelectElement;

    // select.value = null;
    // expect(select.value, null);
    select.value = 'Male';
    expect(select.value, 'Male');
    select.value = 'Female';
    expect(select.value, 'Female');
    select.value = '';
    expect(select.value, '');
  });

  test('fromProperty', async () => {
    const property = DG.Property.js('showPoints', DG.TYPE.STRING, {
      category: 'Fitting',
      description: 'Whether points/candlesticks/none should be rendered',
      defaultValue: 'points', choices: ['points', 'candlesticks', 'both']
    });
    const input = DG.InputBase.forProperty(property, {});
    const view = grok.shell.newView();
    view.root.appendChild(input.root);

    input.value = null;
    expect(input.value == null);
    input.value = 'points';
    expect(input.value, 'points');
    input.value = 'candlesticks';
    expect(input.value, 'candlesticks');
    input.value = 'both';
    expect(input.value, 'both');
    // input.value = '';
    // expect(input.value, '');

    const select = document.querySelector('select.ui-input-editor') as HTMLSelectElement;

    // select.value = null;
    // expect(select.value, null);
    select.value = 'points';
    expect(select.value, 'points');
    select.value = 'candlesticks';
    expect(select.value, 'candlesticks');
    select.value = 'both';
    expect(select.value, 'both');
    select.value = '';
    expect(select.value, '');
  });

  after(async () => {
    grok.shell.closeAll();
  });
}, {owner: 'dkovalyov@datagrok.ai'});

category('UI: Choice input new', () => {
  test('nullable', async () => {
    const t = ui.input.choice('test', { value: '1', items: ['1', '2'], nullable: true });

    const view = grok.shell.newView();
    view.append(t);

    const selector = view.root.querySelector('select.ui-input-editor') as HTMLSelectElement;
    expect(selector.item(0)?.textContent, '');
    expect(selector.item(1)?.textContent, '1');
    expect(selector.item(2)?.textContent, '2');

    expectArray(t.items, [null, '1', '2']);
  });

  test('non-nullable', async () => {
    const t = ui.input.choice('test', { value: '1', items: ['1', '2'], nullable: false });

    const view = grok.shell.newView();
    view.append(t);

    const selector = view.root.querySelector('select.ui-input-editor') as HTMLSelectElement;
    expect(selector.item(0)?.textContent, '1');
    expect(selector.item(1)?.textContent, '2');

    expectArray(t.items, ['1', '2']);
  });

  test('fromFunction', async () => {
    const view = grok.shell.newView();
    const input = ui.input.choice('Sex', { value: 'Male', items: ['Male', 'Female'] });
    view.root.appendChild(input.root);

    input.value = 'Male';
    expect(input.value, 'Male');
    input.value = 'Female';
    expect(input.value, 'Female');
    // input.value = '';
    // expect(input.value, '');

    const select = document.querySelector('select.ui-input-editor') as HTMLSelectElement;

    // select.value = null;
    // expect(select.value, null);
    select.value = 'Male';
    expect(select.value, 'Male');
    select.value = 'Female';
    expect(select.value, 'Female');
    select.value = '';
    expect(select.value, '');
  });

  test('fromProperty', async () => {
    const property = DG.Property.js('showPoints', DG.TYPE.STRING, {
      category: 'Fitting',
      description: 'Whether points/candlesticks/none should be rendered',
      defaultValue: 'points', choices: ['points', 'candlesticks', 'both']
    });
    const input = DG.InputBase.forProperty(property, {});
    const view = grok.shell.newView();
    view.root.appendChild(input.root);

    input.value = null;
    expect(input.value == null);
    input.value = 'points';
    expect(input.value, 'points');
    input.value = 'candlesticks';
    expect(input.value, 'candlesticks');
    input.value = 'both';
    expect(input.value, 'both');
    // input.value = '';
    // expect(input.value, '');

    const select = document.querySelector('select.ui-input-editor') as HTMLSelectElement;

    // select.value = null;
    // expect(select.value, null);
    select.value = 'points';
    expect(select.value, 'points');
    select.value = 'candlesticks';
    expect(select.value, 'candlesticks');
    select.value = 'both';
    expect(select.value, 'both');
    select.value = '';
    expect(select.value, '');
  });

  after(async () => {
    grok.shell.closeAll();
  });
}, {owner: 'dkovalyov@datagrok.ai'});

category('UI: Table input new', () => {
  test('nullable', async () => {
    grok.shell.addTableView(grok.data.demo.demog(10));
    const t = ui.input.table('test', { nullable: true });

    const view = grok.shell.newView();
    view.append(t);

    const selector = view.root.querySelector('select.ui-input-editor') as HTMLSelectElement;
    expect(selector.item(0)?.textContent, '');
    expect(selector.item(1)?.textContent, 'demog 10');
  });

  test('non-nullable', async () => {
    grok.shell.addTableView(grok.data.demo.demog(10));
    const t = ui.input.table('test', { nullable: false });

    const view = grok.shell.newView();
    view.append(t);

    const selector = view.root.querySelector('select.ui-input-editor') as HTMLSelectElement;
    expect(selector.item(0)?.textContent, 'demog 10');
  });

  after(async () => {
    grok.shell.closeAll();
  });
}, {owner: 'dkovalyov@datagrok.ai'});
