import {after, before, category, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {caption, checkHTMLElement, customCaption, stringValue, units} from './utils';

category('UI: Inputs for property', () => {
  const dummyInput = ui.input.int('dummy', {value: 3});
  const props = {
    'stringInput': DG.Property.fromOptions(
      {name: 'stringInput', caption: 'stringCaption', defaultValue: 'defaultString', type: 'string'},
    ),
    'intInput': DG.Property.fromOptions(
      {name: 'intInput', caption: 'intCaption', defaultValue: 3, type: 'int', units: 'kg'},
    ),
    'floatInput': DG.Property.fromOptions(
      {name: 'floatInput', caption: 'doubleCaption', defaultValue: 4.5, type: 'double', units: 'm/s'},
    ),
    'boolInput': DG.Property.fromOptions(
      {name: 'boolInput', caption: 'boolCaption', defaultValue: true, type: 'bool'},
    ),
    'choiceInput': DG.Property.fromOptions(
      {name: 'choiceInput', caption: 'choiceCaption', defaultValue: '1', choices: ['1', '2', '3'], type: 'string'},
    ),
    'tableInput': DG.Property.fromOptions(
      {name: 'tableInput', caption: 'dataframeCaption', defaultValue: null, type: 'dataframe'},
    ),
    // 'columnInput': ui.input.column('', {table: t, value: t.col('age')}),
    // 'columnsInput': ui.input.columns('', {table: t}),
  } as Record<string, DG.Property>;
  let inputs = {} as Record<string, DG.InputBase>;
  let v: DG.View;
  before(async () => {
    inputs = Object.entries(props)
      .map(([inputName, prop]) => ({inputName, input: ui.input.forProperty(prop)}))
      .reduce((acc, config) => {
        acc[config.inputName] = config.input ?? null;

        return acc;
      }, {} as Record<string, DG.InputBase>);
    v = grok.shell.newView('');
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
      customCaption(key, value, v, props[key].caption);
  });

  test('input.stringValue', async () => {
    for (const [key, value] of Object.entries(inputs))
      stringValue(key, value, '.ui-input-editor', v);
  });

  test('input.units.exist', async () => {
    for (const [key, value] of Object.entries(inputs).filter(([key, _]) => props[key].options.units))
      checkHTMLElement(key, value.root, v, '.ui-input-description');
  });

  test('input.units.value', async () => {
    for (const [key, value] of Object.entries(inputs).filter(([key, _]) => props[key].options.units))
      units(key, value, v, '.ui-input-description', props[key].options.units);
  });

  test('input.width', async () => {
    for (const value of Object.values(inputs))
      expect(value.input.clientWidth, dummyInput.input.clientWidth);
  });

  after(async () => {
    grok.shell.closeAll();
  });
});

