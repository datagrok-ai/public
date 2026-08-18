import * as DG from 'datagrok-api/dg';
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes, registerAllFunctions, getRegisteredFuncs} from '../rete/node-factory';
import {propertyChoices, stringChoiceOptions, PropertyPanel} from '../panel/property-panel';
import {getParamDisplayName, getParamDefault, unquoteDefault, getFuncDisplayName} from '../utils/dart-proxy-utils';
import {propertyNameToFriendly} from '../utils/naming';
import {missingRequiredInputs, EXEC_IN_KEY, EXEC_OUT_KEY} from '../rete/scheme';
import {CUSTOM_FUNC_INPUT_EDITORS, CustomInputEditor, customEditorFor} from '../utils/func-input-overrides';
import {EDITOR_SHORTCUT_INPUTS} from '../utils/func-editor-utils';
import {tid} from '../utils/test-ids';
import {makeEditor, destroyEditor, addNode, until} from './test-utils';

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
    expect(stringChoiceOptions(['inner'], false, '')!.join(','), 'inner');
  });

  test('getParamDisplayName falls back to the humanized property name when no caption', async () => {
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
    // Prefer string inputs — no numeric coercion.
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
      const kind = chips!.querySelector('[data-testid="ff-property-type-badge"]') as HTMLElement | null;
      expect(kind?.textContent, 'Function', 'the kind chip leads with a user-facing word');
      const fullName = chips!.querySelector('[data-testid="ff-prop-func-fullname"]') as HTMLElement | null;
      const pkgName = node.dgPackageName ?? '';
      const expectedName = pkgName && node.dgFuncName?.toLowerCase().startsWith(`${pkgName.toLowerCase()}:`) ?
        node.dgFuncName!.slice(pkgName.length + 1) : node.dgFuncName;
      expect(fullName?.textContent, expectedName, 'name chip carries the (deduped) function name');
      if (pkgName) {
        const pkg = chips!.querySelector('[data-testid="ff-prop-func-package"]') as HTMLElement | null;
        expect(pkg?.textContent, pkgName, 'package chip present');
      }
      expect(chips!.querySelectorAll('.funcflow-chip-sep').length >= 1, true,
        'items separated by vertical rules on one line');
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

      const expectedMissing = missingRequiredInputs(join, (k) => e.flow.isInputConnected(join.id, k));
      expect(expectedMissing.length > 0, true, 'JoinTables with one table wired still misses required inputs');
      const rendered = Array.from(panel.root.querySelectorAll('[data-missing]'))
        .map((el) => (el as HTMLElement).dataset.missing);
      expect(rendered.join(','), expectedMissing.join(','), 'missing rows match the helper output in order');

      await e.flow.addConnectionByKeys(tableIn.id, EXEC_OUT_KEY, join.id, EXEC_IN_KEY);
      panel.showNode(join);
      const orderRow = panel.root.querySelector(`[data-conn="${EXEC_IN_KEY}"]`) as HTMLElement | null;
      expect(!!orderRow, true, 'order edge listed in the Run order group');
      expect(orderRow!.textContent!.includes(`runs after ${tableIn.label}`), true,
        `order row names the predecessor (got: "${orderRow!.textContent}")`);

      // Nothing missing → the pane stays collapsed (lazy content), so expand before querying.
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
      expect(node.hiddenInputs.has('subscribeOnChanges'), true, 'registry read onto the node');
      expect(node.hiddenInputs.has('errorBehavior'), true, 'registry read onto the node');
      expect('subscribeOnChanges' in node.inputs, true, 'slot still exists');
      expect('subscribeOnChanges__pt' in node.outputs, true, 'pass-through still exists');
      // The declared default may change — assert the seed follows the declaration without pinning its value.
      expect(typeof node.inputValues['subscribeOnChanges'], 'boolean', 'default still seeded');
      expect(node.passthroughCount, node.dgFunc!.inputs.length, 'pass-through count unchanged');

      // Anchor on the dataframe input — primitive rows are default-hidden.
      await until(() => !!e.container.querySelector(`[data-testid="${tid('socket-input', 'table')}"]`));
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

  test('editor-shortcut inputs render a pencil that opens the function editor', async () => {
    // Driven on `name` via a temporary registration — `expression` has a custom editor, which short-circuits the DG-input path the pencil hangs off.
    const typeName = funcTypeName('AddNewColumn');
    if (!typeName) return;
    const e = makeEditor();
    const panel = new PropertyPanel(e.flow);
    document.body.appendChild(panel.root);
    const shortcut = 'core:addnewcolumn:name';
    EDITOR_SHORTCUT_INPUTS.add(shortcut);
    try {
      const node = await addNode(e.flow, typeName);

      panel.showNode(node);
      const pencilSel = `[data-testid="${tid('prop-input-editor-name')}"]`;
      expect(panel.root.querySelector(pencilSel) == null, true,
        'no pencil when the editor launcher is not wired');

      let openedFor: string | null = null;
      panel.onEditFuncParams = (n): void => {openedFor = n.id;};
      panel.showNode(node);
      const pencil = panel.root.querySelector(pencilSel) as HTMLElement | null;
      expect(!!pencil, true, 'the listed input carries the pencil');
      expect(!!pencil!.closest('[data-param="name"]'), true, 'pencil sits inside its own row');
      pencil!.click();
      expect(openedFor, node.id, 'the pencil opens the editor for the shown node');

      expect(panel.root.querySelector(`[data-testid="${tid('prop-input-editor-type')}"]`) == null, true,
        'off-list inputs get no pencil');
    } finally {
      EDITOR_SHORTCUT_INPUTS.delete(shortcut);
      panel.root.remove();
      destroyEditor(e);
    }
  });

  test('AddNewColumn edits its expression with the inline formula editor', async () => {
    const typeName = funcTypeName('AddNewColumn');
    if (!typeName) return;
    const e = makeEditor();
    const panel = new PropertyPanel(e.flow);
    document.body.appendChild(panel.root);
    try {
      const node = await addNode(e.flow, typeName);
      expect(customEditorFor(node.dgFunc!, 'expression') !== null, true,
        'the expression parameter has a custom editor registered');

      panel.onEditFuncParams = (): void => {};
      node.inputValues['expression'] = '${age} + 1';
      panel.showNode(node);

      const row = panel.root.querySelector('[data-param="expression"]') as HTMLElement | null;
      expect(!!row, true, 'the expression row renders');
      const host = row!.querySelector('.ff-expression-editor') as HTMLElement | null;
      expect(!!host, true, 'and hosts the formula editor');
      expect(host!.getAttribute('data-mode'), 'plain', 'with no table it degrades to a usable input');
      expect((host!.querySelector('input') as HTMLInputElement).value, '${age} + 1',
        'seeded from the stored expression');
      expect(panel.root.querySelector(`[data-testid="${tid('prop-input-editor-expression')}"]`) == null, true,
        'the editor supersedes the pencil');
    } finally {
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

      const factory = customEditorFor(node.dgFunc!, 'fullPath');
      expect(!!factory, true, 'custom editor registered for core:OpenFile fullPath');
      const param = node.dgFunc!.inputs.find((p) => p.name === 'fullPath')!;
      // This editor uses none of the node context, so an inert one is enough.
      const ed = factory!(param, {inputValue: () => undefined, columns: () => null, watch: () => {}});
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

  test('OpenFile: fullPath and sheetName stay off the node body but editable in the panel', async () => {
    const typeName = funcTypeName('OpenFile');
    if (!typeName) return;
    const e = makeEditor();
    const panel = new PropertyPanel(e.flow);
    document.body.appendChild(panel.root);
    try {
      const node = await addNode(e.flow, typeName);
      expect('fullPath' in node.inputs, true, 'fullPath slot still exists');
      expect('sheetName' in node.inputs, true, 'sheetName slot still exists');
      expect(node.hiddenInputs.has('fullPath') && node.hiddenInputs.has('sheetName'), true,
        'both are node-hidden (PANEL_ONLY_FUNC_INPUTS read onto the node)');
      expect(node.requiredInputs.includes('fullPath'), true, 'fullPath still required');

      // Wait for the node card (its real output row) before asserting absence.
      await until(() => !!e.container.querySelector(`[data-testid="${tid('socket-output', 'result')}"]`));
      for (const key of ['fullPath', 'sheetName']) {
        expect(!!e.container.querySelector(`[data-testid="${tid('socket-input', key)}"]`), false,
          `no node row for ${key}`);
        expect(!!e.container.querySelector(`[data-testid="${tid('socket-output', `${key}__pt`)}"]`), false,
          `no pass-through row for ${key}`);
      }

      node.inputValues['fullPath'] = 'System:AppData/Demo/book.xlsx';
      panel.showNode(node);
      expect(!!panel.root.querySelector('[data-param="fullPath"] input'), true, 'panel edits fullPath');
      expect(!!panel.root.querySelector('[data-param="sheetName"] input'), true,
        'panel edits sheetName for an xlsx path');

      node.inputValues['fullPath'] = 'System:AppData/Chem/mol1K.csv';
      panel.showNode(node);
      expect(panel.root.querySelector('[data-param="sheetName"]') == null, true,
        'no Sheet Name row for a non-Excel path');
    } finally {
      panel.root.remove();
      destroyEditor(e);
    }
  });
});
