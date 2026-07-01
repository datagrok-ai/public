/** Inline editing of `string_list` (and its `list<string>` alias) inputs: the
 *  func node seeds the slot editable, the property panel renders a text field,
 *  and the compiler turns the comma-separated value into a JS array of trimmed,
 *  non-empty strings (empty → omitted, so the function's default applies). */
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';
import {registerBuiltinNodes, registerAllFunctions, getRegisteredFuncs} from '../rete/node-factory';
import {emitScript} from '../compiler/script-emitter';
import {PropertyPanel} from '../panel/property-panel';
import {makeEditor, destroyEditor, addNode} from './test-utils';

const SETTINGS = {name: 'T', description: '', tags: []};

/** A registered func with a string_list input, plus that input's name. Aggregate
 *  (fields) is the canonical one; fall back to any func with such an input. */
function stringListFunc(): {typeName: string; param: string} | null {
  const funcs = getRegisteredFuncs();
  const preferred = funcs.find((f) => f.func.name === 'Aggregate' &&
    f.func.inputs.some((i) => String(i.propertyType) === 'string_list'));
  const pick = preferred ?? funcs.find((f) => f.func.inputs.some((i) => String(i.propertyType) === 'string_list'));
  if (!pick) return null;
  const param = pick.func.inputs.find((i) => String(i.propertyType) === 'string_list')!.name;
  return {typeName: pick.nodeTypeName, param};
}

category('Flow: string-list inputs', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('a string_list input is seeded editable (not connection-only)', async () => {
    const f = stringListFunc();
    if (!f) return; // no string_list func on this server — skip
    const e = makeEditor();
    try {
      const node = await addNode(e.flow, f.typeName);
      expect(f.param in node.inputValues, true, `${f.param} is editable`);
      expect(node.inputValues[f.param], '', 'seeded empty');
    } finally {
      destroyEditor(e);
    }
  });

  test('the property panel renders a text field for it', async () => {
    const f = stringListFunc();
    if (!f) return;
    const e = makeEditor();
    const panel = new PropertyPanel(e.flow);
    document.body.appendChild(panel.root);
    try {
      const node = await addNode(e.flow, f.typeName);
      panel.showNode(node);
      const field = panel.root.querySelector(`[data-param="${f.param}"] textarea`);
      expect(!!field, true, 'comma-separated text field present');
    } finally {
      panel.root.remove();
      destroyEditor(e);
    }
  });

  test('a comma-separated value compiles to a trimmed string array', async () => {
    const f = stringListFunc();
    if (!f) return;
    const e = makeEditor();
    try {
      const node = await addNode(e.flow, f.typeName);
      node.inputValues[f.param] = ' alpha, beta ,, gamma ';
      const clean = emitScript(e.flow, SETTINGS);
      expect(clean.includes('["alpha", "beta", "gamma"]'), true, 'trimmed array emitted');
    } finally {
      destroyEditor(e);
    }
  });

  test('an empty string_list value is omitted (function default kept)', async () => {
    const f = stringListFunc();
    if (!f) return;
    const e = makeEditor();
    try {
      const node = await addNode(e.flow, f.typeName);
      node.inputValues[f.param] = '   ,  , ';
      const clean = emitScript(e.flow, SETTINGS);
      expect(new RegExp(`${f.param}\\s*:`).test(clean), false, 'empty list arg not emitted');
    } finally {
      destroyEditor(e);
    }
  });
});
