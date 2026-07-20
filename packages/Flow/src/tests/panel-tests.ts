import * as DG from 'datagrok-api/dg';
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes, registerAllFunctions, getRegisteredFuncs} from '../rete/node-factory';
import {propertyChoices, stringChoiceOptions, PropertyPanel} from '../panel/property-panel';
import {getParamDisplayName, getParamDefault, unquoteDefault, getFuncDisplayName} from '../utils/dart-proxy-utils';
import {propertyNameToFriendly} from '../utils/naming';
import {missingRequiredInputs, EXEC_IN_KEY, EXEC_OUT_KEY} from '../rete/scheme';
import {CUSTOM_FUNC_INPUT_EDITORS, CustomInputEditor, customEditorFor} from '../utils/func-input-overrides';
import {tid} from '../utils/test-ids';
import {makeEditor, destroyEditor, addNode, until} from './test-utils';

/** Registered factory name for a function by name, or null if absent. */
function funcTypeName(name: string): string | null {
  return getRegisteredFuncs().find((f) => f.func.name === name)?.nodeTypeName ?? null;
}

function paneHeaders(root: HTMLElement): string[] {
  return Array.from(root.querySelectorAll('.d4-accordion-pane-header')).map((h) => (h.textContent ?? '').trim());
}

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

  test('getParamDisplayName falls back to the humanized property name when no caption', async () => {
    // No / empty caption → core's propertyNameToFriendly form (the
    // caption-present case is validated end-to-end against live registered
    // funcs in node-factory-tests).
    expect(getParamDisplayName(DG.Property.fromOptions({name: 'minPts', type: 'int'})), 'Min Pts',
      'no caption → the humanized name');
    expect(getParamDisplayName(DG.Property.fromOptions({name: 'minPts', caption: '', type: 'int'})), 'Min Pts',
      'empty caption → the humanized name');
  });

  test('propertyNameToFriendly mirrors core (capitalizeWords ∘ camelCaseToWords ∘ dot→space)', async () => {
    expect(propertyNameToFriendly('maxNumOfSomething'), 'Max Num Of Something');
    expect(propertyNameToFriendly('molBlock'), 'Mol Block');
    expect(propertyNameToFriendly('table'), 'Table');
    expect(propertyNameToFriendly('ratio.split'), 'Ratio Split');
    expect(propertyNameToFriendly('table1'), 'Table1', 'digits do not split (core behavior)');
    expect(propertyNameToFriendly('Log P'), 'Log P', 'already-friendly text passes through');
    expect(propertyNameToFriendly('MW'), 'MW', 'all-caps acronyms preserved (deviation from core)');
    expect(propertyNameToFriendly('HBA'), 'HBA', 'all-caps acronyms preserved');
    expect(propertyNameToFriendly('maxMW'), 'Max MW', 'acronym preserved inside a camelCase name');
    expect(propertyNameToFriendly('RDKitMol'), 'RD Kit Mol', 'acronym-run split keeps the caps');
    expect(propertyNameToFriendly('__exec_in'), '__exec_in', 'non-letter identifiers untouched');
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

  test('header block combines title, chips, and description; func description seeds the input', async () => {
    // Any registered func with its own description exercises the combine.
    const info = getRegisteredFuncs().find((f) => {
      try {
        return !!f.func.description;
      } catch {
        return false;
      }
    });
    if (!info) return; // no described funcs on this stand — skip
    const e = makeEditor();
    const panel = new PropertyPanel(e.flow);
    document.body.appendChild(panel.root);
    try {
      const node = await addNode(e.flow, info.nodeTypeName);
      panel.showNode(node);

      const header = panel.root.querySelector('[data-testid="ff-property-header"]') as HTMLElement | null;
      expect(!!header, true, 'header block rendered');
      expect(!!header!.querySelector('[data-testid="ff-property-title-row"]'), true, 'title row inside the header');
      const descEl = header!.querySelector('[data-param="Description"] textarea') as HTMLTextAreaElement | null;
      expect(!!descEl, true, 'description input inside the header (same padding as the title)');
      expect(descEl!.value, info.func.description, 'empty node description → seeded from the function');

      node.description = 'custom note';
      panel.showNode(node);
      const descEl2 = panel.root.querySelector('[data-param="Description"] textarea') as HTMLTextAreaElement;
      expect(descEl2.value, 'custom note', 'a node-level description wins over the function text');

      expect(paneHeaders(panel.root).includes('Function'), false, 'no separate Function pane anymore');
      const chips = panel.root.querySelector('[data-testid="ff-prop-func-chips"]') as HTMLElement | null;
      expect(!!chips, true, 'chips row rendered for a func node');
      const fullName = chips!.querySelector('[data-testid="ff-prop-func-fullname"]') as HTMLElement | null;
      expect(fullName?.textContent, node.dgFuncName, 'full-name chip carries the qualified name');
      if (node.dgPackageName) {
        const pkg = chips!.querySelector('[data-testid="ff-prop-func-package"]') as HTMLElement | null;
        expect(pkg?.textContent, node.dgPackageName, 'package chip present');
      }
    } finally {
      panel.root.remove();
      destroyEditor(e);
    }
  });

  test('parameters pane is titled with the function display name and expanded', async () => {
    const typeName = funcTypeName('JoinTables');
    if (!typeName) return; // not registered on this server — skip
    const e = makeEditor();
    const panel = new PropertyPanel(e.flow);
    document.body.appendChild(panel.root);
    try {
      const join = await addNode(e.flow, typeName);
      panel.showNode(join);
      const expected = getFuncDisplayName(join.dgFunc!);
      expect(paneHeaders(panel.root).includes(expected), true,
        `pane titled "${expected}" (got: ${paneHeaders(panel.root).join(' | ')})`);
      expect(paneHeaders(panel.root).includes('Input Parameters'), false, 'old pane title gone');
      // Expanded by default → its content (input rows) is already in the DOM.
      expect(!!panel.root.querySelector('[data-param="keys1"]'), true,
        'pane content rendered without a click');
    } finally {
      panel.root.remove();
      destroyEditor(e);
    }
  });

  test('connections pane lists only wired slots and flags missing required ones', async () => {
    const typeName = funcTypeName('JoinTables');
    if (!typeName) return;
    const e = makeEditor();
    const panel = new PropertyPanel(e.flow);
    document.body.appendChild(panel.root);
    try {
      const join = await addNode(e.flow, typeName);
      const tableIn = await addNode(e.flow, 'Inputs/Table Input');
      await e.flow.addConnectionByKeys(tableIn.id, 'table', join.id, 'table1');
      panel.showNode(join);

      const inRow = panel.root.querySelector('[data-conn="table1"]') as HTMLElement | null;
      expect(!!inRow, true, 'wired input listed');
      expect(inRow!.textContent!.includes(`← ${tableIn.label} · Table`), true,
        `input row names its source end, humanized (got: "${inRow!.textContent}")`);
      expect(!!panel.root.querySelector('[data-conn="table2"]'), false, 'unwired input not listed');

      // The Missing group mirrors the shared helper (labels of required inputs
      // neither connected nor filled) — and its presence auto-expands the pane.
      const expectedMissing = missingRequiredInputs(join, (k) => e.flow.isInputConnected(join.id, k));
      expect(expectedMissing.length > 0, true, 'JoinTables with one table wired still misses required inputs');
      const rendered = Array.from(panel.root.querySelectorAll('[data-missing]'))
        .map((el) => (el as HTMLElement).dataset.missing);
      expect(rendered.join(','), expectedMissing.join(','), 'missing rows match the helper output in order');

      // An order edge renders as a plain run-order fact, not a raw __exec row.
      await e.flow.addConnectionByKeys(tableIn.id, EXEC_OUT_KEY, join.id, EXEC_IN_KEY);
      panel.showNode(join);
      const orderRow = panel.root.querySelector(`[data-conn="${EXEC_IN_KEY}"]`) as HTMLElement | null;
      expect(!!orderRow, true, 'order edge listed in the Run order group');
      expect(orderRow!.textContent!.includes(`runs after ${tableIn.label}`), true,
        `order row names the predecessor (got: "${orderRow!.textContent}")`);

      // The wired source node: nothing missing → the pane stays collapsed
      // (content lazily built); expanding it lists the wired output only.
      panel.showNode(tableIn);
      expect(panel.root.querySelectorAll('[data-missing]').length, 0, 'no missing rows on a satisfied node');
      const connHeader = Array.from(panel.root.querySelectorAll('.d4-accordion-pane-header'))
        .find((h) => (h.textContent ?? '').trim() === 'Connections') as HTMLElement | undefined;
      expect(!!connHeader, true, 'Connections pane present');
      connHeader!.click();
      const outRow = panel.root.querySelector('[data-conn="table"]') as HTMLElement | null;
      expect(!!outRow, true, 'source output listed as wired');
      expect(outRow!.textContent!.includes(`→ ${join.label} · `), true,
        `output row names its target end (got: "${outRow!.textContent}")`);
    } finally {
      panel.root.remove();
      destroyEditor(e);
    }
  });

  test('hidden inputs (HIDDEN_FUNC_INPUTS) stay data-carrying but render nowhere', async () => {
    const typeName = funcTypeName('AddNewColumn');
    if (!typeName) return; // not registered on this server — skip
    const e = makeEditor();
    const panel = new PropertyPanel(e.flow);
    document.body.appendChild(panel.root);
    try {
      const node = await addNode(e.flow, typeName);
      // Data layer untouched — slots, seeds, and pass-throughs still exist, so
      // compile and creation-script import/emit round-trip these params.
      expect(node.hiddenInputs.has('subscribeOnChanges'), true, 'registry read onto the node');
      expect(node.hiddenInputs.has('errorBehavior'), true, 'registry read onto the node');
      expect('subscribeOnChanges' in node.inputs, true, 'slot still exists');
      expect('subscribeOnChanges__pt' in node.outputs, true, 'pass-through still exists');
      // GROK-20428 flipped the declared default to true — assert the seed
      // exists and follows the declaration, without pinning its value.
      expect(typeof node.inputValues['subscribeOnChanges'], 'boolean', 'default still seeded');
      expect(node.passthroughCount, node.dgFunc!.inputs.length, 'pass-through count unchanged');

      // Rendered node: rows for hidden inputs (and their pass-throughs) are gone.
      await until(() => !!e.container.querySelector(`[data-testid="${tid('socket-input', 'expression')}"]`));
      expect(!!e.container.querySelector(`[data-testid="${tid('socket-input', 'subscribeOnChanges')}"]`), false,
        'no node row for a hidden input');
      expect(!!e.container.querySelector(`[data-testid="${tid('socket-input', 'errorBehavior')}"]`), false,
        'no node row for a hidden input');
      expect(!!e.container.querySelector(`[data-testid="${tid('socket-output', 'subscribeOnChanges__pt')}"]`), false,
        'no pass-through row either');

      panel.showNode(node);
      expect(!!panel.root.querySelector('[data-param="subscribeOnChanges"]'), false, 'no panel row');
      expect(!!panel.root.querySelector('[data-param="errorBehavior"]'), false, 'no panel row');
      expect(!!panel.root.querySelector('[data-param="expression"]'), true, 'visible param rows still render');
    } finally {
      panel.root.remove();
      destroyEditor(e);
    }
  });

  test('a registered custom editor replaces the default input and routes storage/validity', async () => {
    const typeName = funcTypeName('AddNewColumn');
    if (!typeName) return;
    const e = makeEditor();
    const panel = new PropertyPanel(e.flow);
    document.body.appendChild(panel.root);
    const node = await addNode(e.flow, typeName);
    const nq = node.dgFunc!.nqName;
    const hadEntry = CUSTOM_FUNC_INPUT_EDITORS[nq];
    let valid = true;
    let ed: CustomInputEditor | null = null;
    CUSTOM_FUNC_INPUT_EDITORS[nq] = {
      ...hadEntry,
      name: (): CustomInputEditor => {
        const el = document.createElement('input');
        return ed = {
          element: el,
          getValue: () => el.value,
          setValue: (v) => {el.value = String(v ?? '');},
          isValid: () => valid,
        };
      },
    };
    try {
      node.inputValues['name'] = 'seeded';
      panel.showNode(node);
      const row = panel.root.querySelector('[data-param="name"]') as HTMLElement | null;
      expect(!!row, true, 'row rendered for the overridden input');
      expect(row!.contains(ed!.element), true, 'the row hosts the custom element, not a DG input');
      expect(String(ed!.getValue()), 'seeded', 'editor initialized from inputValues via setValue');

      ed!.onChanged!('x2');
      expect(String(node.inputValues['name']), 'x2', 'a change is stored in inputValues');
      valid = false;
      ed!.onChanged!('bad');
      expect(String(node.inputValues['name']), 'x2', 'isValid() === false blocks the store');
    } finally {
      if (hadEntry) CUSTOM_FUNC_INPUT_EDITORS[nq] = hadEntry;
      else delete CUSTOM_FUNC_INPUT_EDITORS[nq];
      panel.root.remove();
      destroyEditor(e);
    }
  });

  test('OpenFile fullPath renders the file picker; a path round-trips through it', async () => {
    const typeName = funcTypeName('OpenFile');
    if (!typeName) return;
    const e = makeEditor();
    const panel = new PropertyPanel(e.flow);
    document.body.appendChild(panel.root);
    try {
      const node = await addNode(e.flow, typeName);
      panel.showNode(node);
      const row = panel.root.querySelector('[data-param="fullPath"]') as HTMLElement | null;
      expect(!!row, true, 'fullPath row rendered');
      expect(!!row!.querySelector('input'), true, 'file editor input inside the row');

      // The registered editor against the live server: setValue resolves the
      // path (existence-checked), onChanged reports the plain full-path string.
      const factory = customEditorFor(node.dgFunc!, 'fullPath');
      expect(!!factory, true, 'custom editor registered for core:OpenFile fullPath');
      const param = node.dgFunc!.inputs.find((p) => p.name === 'fullPath')!;
      const ed = factory!(param);
      let reported: unknown = null;
      ed.onChanged = (v): void => {reported = v;};
      const path = 'System:AppData/Chem/mol1K.csv';
      ed.setValue(path);
      await until(() => reported !== null, 15000);
      expect(String(reported), path, 'onChanged reports the full-path string');
      expect(String(ed.getValue()), path, 'getValue returns the resolved full path');
    } finally {
      panel.root.remove();
      destroyEditor(e);
    }
  });
});
