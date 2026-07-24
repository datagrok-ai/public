/** Configured input-node values: resolution, header emission, the autorun
 *  gate they unblock, the on-node editor, and the blocked-autorun report. */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes, createNode} from '../rete/node-factory';
import {resolveInputValue, inputBlockReason} from '../utils/input-values';
import {emitScript} from '../compiler/script-emitter';
import {ExecutionController} from '../execution/execution-controller';
import {NodeExecStatus} from '../execution/execution-state';
import {InputValueControl} from '../rete/nodes/input-value-control';
import {FlowEditor, GraphEdit} from '../rete/flow-editor';
import {FlowNode} from '../rete/scheme';
import {makeEditor, destroyEditor, addNode, until} from './test-utils';

const SETTINGS = {name: 'InputValueTest', description: '', tags: []};

function make(typeName: string): FlowNode {
  registerBuiltinNodes();
  const n = createNode(typeName);
  if (!n) throw new Error(`Unknown node type: ${typeName}`);
  return n;
}

category('Flow: input values', () => {
  test('resolveInputValue: scalars, lists, map', async () => {
    const int = make('Inputs/Int Input');
    expect(resolveInputValue(int).ok, true, 'int defaults to 0 — configured');
    expect(resolveInputValue(int).value, 0, 'int zero default');
    int.properties['defaultValue'] = '7.6';
    expect(resolveInputValue(int).value, 8, 'int rounds');
    int.properties['defaultValue'] = 'abc';
    expect(resolveInputValue(int).ok, false, 'non-numeric int is unresolvable');

    const bool = make('Inputs/Boolean Input');
    expect(resolveInputValue(bool).ok, true, 'bool always resolvable');
    expect(resolveInputValue(bool).value, false, 'bool default false');

    const str = make('Inputs/String Input');
    expect(resolveInputValue(str).ok, false, 'empty string blocks (nothing chosen yet)');
    str.properties['nullable'] = true;
    expect(resolveInputValue(str).ok, true, 'nullable empty string is a value');
    str.properties['nullable'] = false;
    str.properties['defaultValue'] = 'hello';
    expect(resolveInputValue(str).value, 'hello', 'string value');

    const dt = make('Inputs/DateTime Input');
    dt.properties['defaultValue'] = 'not a date';
    expect(resolveInputValue(dt).ok, false, 'garbage date is unresolvable');
    dt.properties['defaultValue'] = '2026-07-23T10:00:00.000Z';
    expect(resolveInputValue(dt).ok, true, 'ISO date resolves');

    const list = make('Inputs/String List Input');
    list.properties['defaultValue'] = 'a, b , ,c';
    expect(JSON.stringify(resolveInputValue(list).value), '["a","b","c"]', 'comma list trimmed');

    const map = make('Inputs/Map Input');
    map.properties['defaultValue'] = '{"k": 1}';
    expect((resolveInputValue(map).value as {k: number}).k, 1, 'map parses JSON');
    map.properties['defaultValue'] = '{nope';
    expect(resolveInputValue(map).ok, false, 'broken JSON is unresolvable');

    const blob = make('Inputs/Blob Input');
    expect(resolveInputValue(blob).ok, false, 'blob has no inline value');
  });

  test('resolveInputValue: dataframe by name and by live reference', async () => {
    const node = make('Inputs/Table Input');
    expect(resolveInputValue(node).ok, false, 'no table selected');
    node.properties['defaultValue'] = 'ff test table';
    expect(resolveInputValue(node).ok, false, 'named table is not open');
    expect(resolveInputValue(node).reason?.includes('ff test table'), true, 'reason names the table');

    const df = DG.DataFrame.fromCsv('x\n1\n2');
    df.name = 'ff test table';
    grok.shell.addTable(df);
    try {
      const byName = resolveInputValue(node);
      expect(byName.ok, true, 'open table resolves by name');
      expect((byName.value as DG.DataFrame).rowCount, 2, 'the right table');

      const other = DG.DataFrame.fromCsv('y\n1');
      node.transientValue = other; // an uploaded-into-the-input table
      expect(resolveInputValue(node).value === other, true, 'live reference wins over the name');
    } finally {
      node.transientValue = undefined;
      grok.shell.closeTable(df);
    }
  });

  test('inputBlockReason names the node and the parameter', async () => {
    const node = make('Inputs/Table Input');
    const reason = inputBlockReason(node)!;
    expect(reason.includes('Table Input'), true, 'names the node');
    expect(reason.includes('df'), true, 'names the param');
    const int = make('Inputs/Int Input');
    expect(inputBlockReason(int), null, 'a configured input does not block');
  });

  test('header emission: scalar defaults stay, table names stay out', async () => {
    const e = makeEditor();
    try {
      const table = await addNode(e.flow, 'Inputs/Table Input');
      table.properties['defaultValue'] = 'demog';
      const int = await addNode(e.flow, 'Inputs/Int Input');
      int.properties['defaultValue'] = 5;
      const out = await addNode(e.flow, 'Outputs/Table Output');
      await e.flow.addConnectionByKeys(table.id, 'table', out.id, 'table');
      const script = emitScript(e.flow, SETTINGS, {});
      expect(script.includes('//input: dataframe df\n'), true, 'no default on the dataframe header');
      expect(script.includes('//input: dataframe df ='), false, 'table name never leaks into the header');
      expect(script.includes('//input: int n = 5'), true, 'scalar default still emitted');
    } finally {
      destroyEditor(e);
    }
  });

  test('runAutorun: unconfigured table input skips, configured one runs silently', async () => {
    const e = makeEditor();
    const df = DG.DataFrame.fromCsv('a\n1\n2\n3\n');
    df.name = 'ff autorun table';
    try {
      const input = await addNode(e.flow, 'Inputs/Table Input');
      const output = await addNode(e.flow, 'Outputs/Table Output');
      await e.flow.addConnectionByKeys(input.id, 'table', output.id, 'table');
      const ctrl = new ExecutionController(e.flow);
      expect(ctrl.runAutorun(new Set([output.id]), SETTINGS), 'skipped',
        'an UNCONFIGURED input in the run set would open a dialog — autorun must not');

      grok.shell.addTable(df);
      input.properties['defaultValue'] = df.name;
      expect(ctrl.runAutorun(new Set([output.id]), SETTINGS), 'started',
        'a configured input feeds the prepared call — autorun proceeds');
      const done = await until(() =>
        ctrl.state.getNodeState(output.id)?.status === NodeExecStatus.completed, 8000);
      expect(done, true, 'the run completed');
      expect(document.querySelector('.d4-dialog') == null, true, 'no parameter dialog appeared');
    } finally {
      grok.shell.closeTable(df);
      destroyEditor(e);
    }
  });

  test('explicit Run with configured scalar input needs no dialog', async () => {
    const e = makeEditor();
    try {
      const input = await addNode(e.flow, 'Inputs/Int Input');
      input.properties['defaultValue'] = 42;
      const output = await addNode(e.flow, 'Outputs/Value Output');
      await e.flow.addConnectionByKeys(input.id, 'value', output.id, 'value');
      const ctrl = new ExecutionController(e.flow);
      ctrl.runInstrumented(SETTINGS);
      const done = await until(() =>
        ctrl.state.getNodeState(output.id)?.status === NodeExecStatus.completed, 8000);
      expect(done, true, 'the run completed without asking');
      expect(document.querySelector('.d4-dialog') == null, true, 'no parameter dialog appeared');
    } finally {
      destroyEditor(e);
    }
  });

  test('autorunBlockers: lists unconfigured pending inputs, empty when configured', async () => {
    const e = makeEditor();
    try {
      const input = await addNode(e.flow, 'Inputs/Table Input');
      const output = await addNode(e.flow, 'Outputs/Table Output');
      await e.flow.addConnectionByKeys(input.id, 'table', output.id, 'table');
      const ctrl = new ExecutionController(e.flow);
      const blockers = ctrl.autorunBlockers();
      expect(blockers.length, 1, 'one blocker');
      expect(blockers[0].includes('Table Input'), true, 'the blocker names the node');

      const int = make('Inputs/Int Input'); // configured by default
      await e.flow.addNodeAt(int, 300, 0);
      expect(ctrl.autorunBlockers().length, 1, 'configured inputs never block');
    } finally {
      destroyEditor(e);
    }
  });

  test('node body renders the DG value editor; panel sync round-trips', async () => {
    const e = makeEditor();
    try {
      const int = await addNode(e.flow, 'Inputs/Int Input');
      const rendered = await until(() =>
        e.container.querySelector('[data-testid="ff-node-value-input"] .ui-input-root') != null, 4000);
      expect(rendered, true, 'the node body hosts a DG input');

      // A store edit made elsewhere (the context panel path goes through
      // notifyNodeParamsChanged) must show up in the on-node editor.
      int.properties['defaultValue'] = 123;
      e.flow.notifyNodeParamsChanged(int.id);
      const ctl = int.controls['value'] as InputValueControl;
      expect(ctl instanceof InputValueControl, true, 'the value control is registered');
      const host = e.container.querySelector('[data-testid="ff-node-value-input"]') as HTMLElement;
      const synced = await until(() => {
        const editor = host.querySelector('input') as HTMLInputElement | null;
        return editor?.value === '123';
      }, 3000);
      expect(synced, true, 'panel-side edits sync into the node editor');
    } finally {
      destroyEditor(e);
    }
  });

  test('editing through the editor reports params-changed exactly once', async () => {
    registerBuiltinNodes();
    const edits: GraphEdit[] = [];
    const container = document.createElement('div');
    container.style.cssText = 'width:1000px;height:700px;position:absolute;left:-10000px;';
    document.body.appendChild(container);
    const flow = new FlowEditor(container, {onGraphEdited: (edit) => edits.push(edit)});
    try {
      const int = createNode('Inputs/Int Input')!;
      await flow.addNodeAt(int, 0, 0);
      await until(() => container.querySelector('[data-testid="ff-node-value-input"] input') != null, 4000);
      const before = edits.filter((x) => x.kind === 'params-changed').length;
      expect(before, 0, 'rendering the editor is not an edit');

      const ctl = int.controls['value'] as InputValueControl;
      // A programmatic sync (same value) must not report either.
      ctl.sync();
      expect(edits.filter((x) => x.kind === 'params-changed').length, 0, 'sync is not an edit');

      // A real user edit through the DG input reports once.
      const editor = container.querySelector('[data-testid="ff-node-value-input"] input') as HTMLInputElement;
      editor.value = '55';
      editor.dispatchEvent(new Event('input', {bubbles: true}));
      editor.dispatchEvent(new Event('change', {bubbles: true}));
      const reported = await until(() =>
        edits.filter((x) => x.kind === 'params-changed').length === 1, 3000);
      expect(reported, true, 'exactly one params-changed for one edit');
      expect(String(int.properties['defaultValue']), '55', 'the store took the value');
    } finally {
      flow.destroy();
      container.remove();
    }
  });
});
