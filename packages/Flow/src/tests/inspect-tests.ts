/** Tests for the "inspect anywhere" keystone (slice-compile) and the pre-run
 *  "Needs input" hint helper. */
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes, createNode} from '../rete/node-factory';
import {FlowNode, missingRequiredInputs} from '../rete/scheme';
import {sliceUpTo} from '../compiler/graph-compiler';
import {emitScript} from '../compiler/script-emitter';
import {makeEditor, destroyEditor, addNode, TestEditor} from './test-utils';

category('Flow: inspect / slice', () => {
  before(async () => {
    registerBuiltinNodes();
  });

  test('sliceUpTo gathers a node and all its upstream ancestors', async () => {
    const e: TestEditor = makeEditor();
    try {
      // const(a) -> ToString(b) -> ValueOutput(out);  separate const(c) (not upstream of b)
      const a = await addNode(e.flow, 'Constants/String');
      const b = await addNode(e.flow, 'Utilities/ToString');
      const out = await addNode(e.flow, 'Outputs/Value Output');
      const c = await addNode(e.flow, 'Constants/Int');
      await e.flow.addConnectionByKeys(a.id, 'value', b.id, 'value');
      await e.flow.addConnectionByKeys(b.id, 'text', out.id, 'value');

      const sliceB = sliceUpTo(e.flow, b.id);
      expect(sliceB.has(b.id), true, 'slice includes the target');
      expect(sliceB.has(a.id), true, 'slice includes the ancestor');
      expect(sliceB.has(out.id), false, 'slice excludes a downstream node');
      expect(sliceB.has(c.id), false, 'slice excludes an unrelated node');

      const sliceOut = sliceUpTo(e.flow, out.id);
      expect(sliceOut.has(a.id) && sliceOut.has(b.id) && sliceOut.has(out.id), true, 'full chain to the output');
      expect(sliceOut.has(c.id), false, 'still excludes the unrelated const');
    } finally {
      destroyEditor(e);
    }
  });

  test('emitScript with onlyNodeIds drops out-of-slice steps', async () => {
    const e: TestEditor = makeEditor();
    try {
      const a = await addNode(e.flow, 'Constants/String');
      const b = await addNode(e.flow, 'Utilities/ToString');
      const out = await addNode(e.flow, 'Outputs/Value Output');
      await e.flow.addConnectionByKeys(a.id, 'value', b.id, 'value');
      await e.flow.addConnectionByKeys(b.id, 'text', out.id, 'value');

      const slice = sliceUpTo(e.flow, b.id);
      const settings = {name: 'S', description: '', tags: []};
      const sliced = emitScript(e.flow, settings, {onlyNodeIds: slice});
      const full = emitScript(e.flow, settings);

      // The sliced script stops at ToString; the full one declares the output param.
      expect(/\/\/output:/.test(full), true, 'full script declares an output');
      expect(/\/\/output:/.test(sliced), false, 'sliced script has no output node');
    } finally {
      destroyEditor(e);
    }
  });

  test('missingRequiredInputs reports unfilled, unconnected required inputs', async () => {
    // A bare FlowNode with two declared required inputs; neither connected.
    const node = new FlowNode('probe');
    node.requiredInputs = ['table', 'molecules'];
    (node.inputs as Record<string, {label?: string}>)['table'] = {label: 'table'};
    (node.inputs as Record<string, {label?: string}>)['molecules'] = {label: 'molecules'};

    // Nothing connected, nothing filled → both missing.
    let missing = missingRequiredInputs(node, () => false);
    expect(missing.length, 2);

    // Fill one via an input value → only the other remains.
    node.inputValues['molecules'] = 'smiles';
    missing = missingRequiredInputs(node, () => false);
    expect(missing.join(','), 'table');

    // Connect the remaining one → none missing.
    missing = missingRequiredInputs(node, (k) => k === 'table');
    expect(missing.length, 0);
  });
});
