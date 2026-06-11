import * as DG from 'datagrok-api/dg';
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {
  registerBuiltinNodes, registerAllFunctions, createNode, ensureFuncNodeType, getRegisteredTypeNames,
} from '../rete/node-factory';
import {FuncNode} from '../rete/nodes/func-node';

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
});
