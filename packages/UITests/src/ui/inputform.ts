import {after, assure, before, category, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {checkHTMLElement, customCaption, stringValue, units} from './utils';

category('UI: Inputs via InputForm', () => {
  const dummyInput = ui.input.int('dummy', {value: 3});
  let inputs = {} as Record<string, DG.InputBase>;
  let funcCall: DG.FuncCall;
  let form: DG.InputForm;
  let v: DG.View;
  before(async () => {
    funcCall = (await grok.functions.eval('UITests:Dummy')).prepare({
      'stringInput': 'stringDefault',
      'intInput': 3,
      'doubleInput': 4.5,
      'boolInput': true,
      'choiceInput': '1',
      'tableInput': grok.data.demo.demog(10),
    });

    form = await DG.InputForm.forFuncCall(funcCall);
    inputs = {
      'stringInput': form.getInput('stringInput') ?? null,
      'intInput': form.getInput('intInput') ?? null,
      'doubleInput': form.getInput('doubleInput') ?? null,
      'boolInput': form.getInput('boolInput') ?? null,
      'choiceInput': form.getInput('choiceInput') ?? null,
      'tableInput': form.getInput('tableInput') ?? null,
      // 'columnInput': ui.input.column('', {table: t, value: t.col('age')}),
      // 'columnsInput': ui.input.columns('', {table: t}),
    };
    v = grok.shell.newView('');
  });

  test('inputform', async () => {
    assure.notNull(form);
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

  test('input.customCaption', async () => {
    for (const [key, value] of Object.entries(inputs))
      customCaption(key, value, v, funcCall.inputParams[key].property.caption, true);
  });

  test('input.stringValue', async () => {
    for (const [key, value] of Object.entries(inputs))
      stringValue(key, value, '.ui-input-editor', v);
  });

  test('input.units.exist', async () => {
    for (const [key, value] of Object.entries(inputs)
      .filter(([key, _]) => ['intInput', 'doubleInput'].includes(key)))
      checkHTMLElement(key, value.root, v, '.ui-input-description');
  });

  test('input.units.value', async () => {
    for (const [key, value] of Object.entries(inputs)
      .filter(([key, _]) => ['intInput', 'doubleInput'].includes(key)))
      units(key, value, v, '.ui-input-description', funcCall.inputParams[key].property.options.units);
  });

  test('input.width', async () => {
    for (const value of Object.values(inputs))
      expect(value.input.clientWidth, dummyInput.input.clientWidth);
  });

  after(async () => {
    grok.shell.closeAll();
  });
}, {owner: 'dkovalyov@datagrok.ai'});
