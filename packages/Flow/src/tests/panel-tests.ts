import * as DG from 'datagrok-api/dg';
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes, registerAllFunctions, getRegisteredFuncs} from '../rete/node-factory';
import {propertyChoices, stringChoiceOptions, PropertyPanel} from '../panel/property-panel';
import {getParamDisplayName, getParamDefault, unquoteDefault} from '../utils/dart-proxy-utils';
import {makeEditor, destroyEditor, addNode} from './test-utils';

category('Flow: property panel', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('stringChoiceOptions: no choices → null (free-text field)', async () => {
    expect(stringChoiceOptions([], false, ''), null);
    expect(stringChoiceOptions([], true, 'x'), null);
  });

  test('stringChoiceOptions: choices kept as-is when not nullable', async () => {
    const opts = stringChoiceOptions(['inner', 'outer'], false, 'inner');
    expect(opts!.join(','), 'inner,outer');
  });

  test('stringChoiceOptions: nullable prepends an empty option', async () => {
    const opts = stringChoiceOptions(['inner', 'outer'], true, 'inner');
    expect(opts!.join(','), ',inner,outer');
  });

  test('stringChoiceOptions: a current value outside the choices is preserved', async () => {
    const opts = stringChoiceOptions(['inner', 'outer'], false, 'custom');
    expect(opts!.join(','), 'custom,inner,outer');
    // Empty current is not added as a bogus option.
    expect(stringChoiceOptions(['inner'], false, '')!.join(','), 'inner');
  });

  test('getParamDisplayName falls back to the property name when no caption', async () => {
    // No / empty caption → the identity name (the caption-present case is
    // validated end-to-end against live registered funcs in node-factory-tests).
    expect(getParamDisplayName(DG.Property.fromOptions({name: 'minPts', type: 'int'})), 'minPts',
      'no caption → the property name');
    expect(getParamDisplayName(DG.Property.fromOptions({name: 'minPts', caption: '', type: 'int'})), 'minPts',
      'empty caption → the property name');
  });

  test('unquoteDefault strips one pair of wrapping quotes', async () => {
    expect(unquoteDefault(`'inner'`), 'inner');
    expect(unquoteDefault('"inner"'), 'inner');
    expect(unquoteDefault('plain'), 'plain');
    expect(unquoteDefault(`  'padded'  `), 'padded');
    expect(unquoteDefault(`'`), `'`); // a lone quote is not a pair
    expect(unquoteDefault(''), '');
  });

  test('a declared default seeds the node and initializes the panel editor', async () => {
    // Find a registered func with a primitive input that declares a default
    // (defaultValue ?? initialValue) — prefer strings (no numeric coercion).
    let found: {typeName: string; param: string; def: unknown} | null = null;
    for (const pass of [['string'], ['int', 'double', 'num']]) {
      for (const info of getRegisteredFuncs()) {
        for (const p of info.func.inputs) {
          if (!pass.includes(String(p.propertyType))) continue;
          const def = getParamDefault(p);
          if (def === undefined || String(def) === '') continue;
          found = {typeName: info.nodeTypeName, param: p.name, def};
          break;
        }
        if (found) break;
      }
      if (found) break;
    }
    if (!found) return; // no declared defaults on this stand — skip

    const e = makeEditor();
    const panel = new PropertyPanel(e.flow);
    document.body.appendChild(panel.root);
    try {
      const node = await addNode(e.flow, found.typeName);
      const seeded = node.inputValues[found.param];
      expect(String(seeded), String(found.def), `inputValues seeded with the default (${found.typeName})`);

      panel.showNode(node);
      const sel = `[data-param="${found.param}"] input, [data-param="${found.param}"] select, ` +
        `[data-param="${found.param}"] textarea`;
      const editor = panel.root.querySelector(sel) as HTMLInputElement | null;
      expect(!!editor, true, 'editor rendered');
      expect(editor!.value, String(seeded), 'editor initialized with the seeded default');
    } finally {
      panel.root.remove();
      destroyEditor(e);
    }
  });

  test('propertyChoices reads declared choices from a live function input', async () => {
    // Find any registered function input that declares choices (e.g. JoinTables
    // joinType). Validates the Dart-proxy `choices` read end to end.
    let found: {func: string; param: string; choices: string[]} | null = null;
    for (const info of getRegisteredFuncs()) {
      for (const inp of info.func.inputs) {
        if (String(inp.propertyType) !== 'string') continue;
        const choices = propertyChoices(inp);
        if (choices.length > 0) {found = {func: info.func.name, param: inp.name, choices}; break;}
      }
      if (found) break;
    }
    if (!found) return; // no choice-bearing string inputs on this server — skip
    expect(found.choices.length > 0, true, `expected choices for ${found.func}.${found.param}`);
    expect(found.choices.every((c) => typeof c === 'string' && c.length > 0), true, 'choices are non-empty strings');
  });
});
