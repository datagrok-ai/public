import * as DG from 'datagrok-api/dg';
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {
  registerBuiltinNodes, registerAllFunctions, createNode, ensureFuncNodeType, getRegisteredTypeNames,
  getRegisteredFuncs,
} from '../rete/node-factory';
import {FuncNode} from '../rete/nodes/func-node';
import {getParamDescription} from '../utils/dart-proxy-utils';

category('Flow: node-factory', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('createNode builds known built-in types', async () => {
    const table = createNode('Inputs/Table Input');
    expect(table !== null, true);
    expect(table!.dgNodeType, 'input');
    expect(table!.dgTypeName, 'Inputs/Table Input'); // stamped for serialization
    expect('table' in table!.outputs, true);

    const out = createNode('Outputs/Table Output');
    expect(out!.dgNodeType, 'output');
    expect('table' in out!.inputs, true);

    const str = createNode('Constants/String');
    expect('value' in str!.outputs, true);
  });

  test('createNode returns null for unknown type', async () => {
    expect(createNode('Nope/Does Not Exist'), null);
  });

  test('built-in registry includes the expected core types', async () => {
    const names = new Set(getRegisteredTypeNames());
    for (const t of ['Inputs/Table Input', 'Outputs/Table Output', 'Outputs/Value Output',
      'Constants/String', 'Constants/Int', 'Constants/Boolean'])
      expect(names.has(t), true, `registry has ${t}`);
  });

  test('ensureFuncNodeType is idempotent for the same function', async () => {
    const func = DG.Func.find({name: 'AddNewColumn'})[0] ?? DG.Func.find({})[0];
    expect(func != null, true);
    const a = ensureFuncNodeType(func);
    const b = ensureFuncNodeType(func);
    expect(a, b); // same qualified name → same node type, no duplicate registry entry
    const node = createNode(a);
    expect(node instanceof FuncNode, true);
    expect(node!.dgFunc?.name, func.name);
  });

  test('FuncNode builds a pass-through output per input', async () => {
    const func = DG.Func.find({name: 'AddNewColumn'})[0];
    if (!func) return; // AddNewColumn should always exist, but guard anyway
    const node = new FuncNode(func);
    expect(node.passthroughCount, func.inputs.length);
    for (const inp of func.inputs)
      expect(`${inp.name}__pt` in node.outputs, true, `pass-through for ${inp.name}`);
  });

  test('FuncNode records the source package (shown in the context panel)', async () => {
    const info = getRegisteredFuncs().find((f) => f.packageName === 'Chem') ??
      getRegisteredFuncs().find((f) => !!f.packageName);
    if (!info) return; // no package funcs on this stand — skip
    const node = new FuncNode(info.func);
    expect(node.dgPackageName, info.packageName, 'package name captured on the node');
  });

  test('FuncNode captures param descriptions for socket/panel tooltips', async () => {
    // Find any registered func whose inputs OR outputs declare a description
    // (via @grok.decorators.param) — validates the Dart-proxy read end to end.
    let found: {func: DG.Func; param: string} | null = null;
    for (const info of getRegisteredFuncs()) {
      for (const p of [...info.func.inputs, ...info.func.outputs])
        if (getParamDescription(p)) {found = {func: info.func, param: p.name}; break;}
      if (found) break;
    }
    if (!found) return; // no described params on this stand — skip
    const node = new FuncNode(found.func);
    const all = {...node.inputDescriptions, ...node.outputDescriptions};
    expect(found.param in all, true, `description captured for ${found.func.name}.${found.param}`);
    expect(all[found.param].length > 0, true, 'description is non-empty');
  });

  test('per-function color override pins SetVar to red', async () => {
    const setVar = DG.Func.find({name: 'SetVar'})[0];
    if (!setVar) return;
    const node = new FuncNode(setVar);
    expect((node as unknown as {color?: string}).color, '#EF5350', 'SetVar title is red');

    // A function without an override falls back to role/default coloring.
    const other = DG.Func.find({name: 'AddNewColumn'})[0];
    if (other)
      expect((new FuncNode(other) as unknown as {color?: string}).color !== '#EF5350', true);
  });
});
