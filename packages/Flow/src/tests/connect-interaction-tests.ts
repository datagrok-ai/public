/** Drop-on-node connection shortcut: dropping an output drag anywhere on a node
 *  with a single compatible, unwired input connects to it. `soleCompatibleInput`
 *  is the decision behind that — exercised here on the data layer (no DOM). */
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';
import {registerBuiltinNodes, registerAllFunctions, getRegisteredFuncs} from '../rete/node-factory';
import {makeEditor, destroyEditor, addNode} from './test-utils';

function typeNameOf(funcName: string): string | null {
  return getRegisteredFuncs().find((f) => f.func.name === funcName)?.nodeTypeName ?? null;
}

category('Flow: connect interaction', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('one compatible free input → that key (drop-on-node connects)', async () => {
    const addCol = typeNameOf('AddNewColumn');
    if (!addCol) return; // not on this server — skip
    const e = makeEditor();
    try {
      const src = await addNode(e.flow, 'Inputs/Table Input');   // output 'table' : dataframe
      const node = await addNode(e.flow, addCol);                // only dataframe input is 'table'
      expect(e.flow.soleCompatibleInput(src.id, 'table', node.id), 'table', 'sole input resolved');
    } finally {
      destroyEditor(e);
    }
  });

  test('several compatible free inputs → null (must aim at a pin)', async () => {
    const join = typeNameOf('JoinTables');
    if (!join) return;
    const e = makeEditor();
    try {
      const src = await addNode(e.flow, 'Inputs/Table Input');
      const node = await addNode(e.flow, join);                  // table1 AND table2 are dataframe
      expect(e.flow.soleCompatibleInput(src.id, 'table', node.id), null, 'ambiguous → null');
      // Wiring table1 leaves table2 as the only free compatible input.
      await e.flow.addConnectionByKeys(src.id, 'table', node.id, 'table1');
      expect(e.flow.soleCompatibleInput(src.id, 'table', node.id), 'table2', 'remaining free input');
    } finally {
      destroyEditor(e);
    }
  });

  test('no compatible input → null (nothing to connect)', async () => {
    const e = makeEditor();
    try {
      const src = await addNode(e.flow, 'Inputs/Table Input');   // dataframe output
      const target = await addNode(e.flow, 'Inputs/String Input'); // an input node — no data inputs
      expect(e.flow.soleCompatibleInput(src.id, 'table', target.id), null, 'no inputs → null');
    } finally {
      destroyEditor(e);
    }
  });

  test('the source node is never its own target', async () => {
    const join = typeNameOf('JoinTables');
    if (!join) return;
    const e = makeEditor();
    try {
      const node = await addNode(e.flow, join);
      // Same node id as source and target — even though it has free dataframe
      // inputs, a node can't connect to itself (caller guards on id, but the
      // query is harmless either way).
      const key = e.flow.soleCompatibleInput(node.id, 'result', node.id);
      // 'result' (dataframe) would match table1/table2, but two free inputs → null.
      expect(key, null);
    } finally {
      destroyEditor(e);
    }
  });
});
