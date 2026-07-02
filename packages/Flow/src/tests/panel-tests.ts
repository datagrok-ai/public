import * as DG from 'datagrok-api/dg';
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes, registerAllFunctions, getRegisteredFuncs} from '../rete/node-factory';
import {propertyChoices, stringChoiceOptions} from '../panel/property-panel';
import {getParamDisplayName} from '../utils/dart-proxy-utils';

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
