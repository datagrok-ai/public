/** Inputs annotated `{editor: Package:func}`: the editor function's input renders in the node's panel; a flow
 *  input wired into such an input is edited there and keeps the value; the run dialog keeps the plain input. */

import * as DG from 'datagrok-api/dg';
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes, registerAllFunctions, ensureFuncNodeType} from '../rete/node-factory';
import {customEditorFor, CustomEditorContext} from '../utils/func-input-overrides';
import {editorAnnotation} from '../panel/editors/annotated-editor';
import {PropertyPanel} from '../panel/property-panel';
import {emitScript} from '../compiler/script-emitter';
import {resolveInputValue} from '../utils/input-values';
import {makeEditor, destroyEditor, addNode, until} from './test-utils';
import {FlowEditor} from '../rete/flow-editor';
import {InputValueControl} from '../rete/nodes/input-value-control';
import * as ui from 'datagrok-api/ui';

const SETTINGS = {name: 'AnnotatedEditorTest', description: '', tags: []};
/** A registered function returning a `DG.InputBase` (PowerPack's cron value editor); the tests needing it skip
 *  without it. */
const EDITOR = 'PowerPack:cronInput';

function editorDeployed(): boolean {
  return DG.Func.find({package: 'PowerPack', name: 'cronInput'}).length > 0;
}

function annotatedScript(editorRef: string): DG.Script {
  return DG.Script.create('//name: ffAnnotatedEditorTest\n//language: javascript\n//input: dataframe table\n' +
    `//input: string schedule {editor: ${editorRef}}\n//output: string out\nout = schedule;\n`);
}

function context(overrides: Partial<CustomEditorContext> = {}): CustomEditorContext {
  return {inputValue: () => undefined, columns: () => null, table: () => null, isConnected: () => true,
    watch: () => {}, ...overrides};
}

category('Flow: editor annotation', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('an input naming an editor gets it; the others keep theirs', async () => {
    const script = annotatedScript(EDITOR);
    const schedule = script.inputs.find((p) => p.name === 'schedule')!;
    expect(editorAnnotation(schedule), EDITOR, 'the annotation is read off the parameter');
    expect(customEditorFor(script, 'schedule') != null, true, 'the annotated input has an editor');
    expect(customEditorFor(script, 'table'), null, 'an input without one has none');
  });

  test('the editor is what the function returns, and its edits come back as the stored value', async () => {
    if (!editorDeployed()) return;
    const script = annotatedScript(EDITOR);
    const param = script.inputs.find((p) => p.name === 'schedule')!;
    const ed = customEditorFor(script, 'schedule')!(param, context());
    document.body.appendChild(ed.element);
    try {
      const changes: unknown[] = [];
      ed.onChanged = (v): void => {changes.push(v);};
      ed.setValue('0 * * * *');
      const select = (): HTMLSelectElement | null => ed.element.querySelector('select.pp-cron-preset-select');
      expect(await until(() => select() != null, 10_000), true, 'the function\'s input is mounted');
      expect(ed.getValue(), '0 * * * *', 'it shows the stored value');
      select()!.value = '*/5 * * * *';
      select()!.dispatchEvent(new Event('change', {bubbles: true}));
      expect(await until(() => changes.includes('*/5 * * * *')), true, 'an edit reports the new value');
    } finally {
      ed.detach?.();
      ed.element.remove();
    }
  });

  test('an editor that is not registered, or whose tables are not there yet, says so', async () => {
    const missing = annotatedScript('NoSuchPackage:noSuchEditor');
    const param = missing.inputs.find((p) => p.name === 'schedule')!;
    const none = customEditorFor(missing, 'schedule')!(param, context());
    expect(none.element.textContent?.includes('not registered'), true, 'an unknown editor is named');

    // AddNewColumn takes a table: before the flow has run, the panel offers to run up to the node.
    const tableEditor = annotatedScript('AddNewColumn');
    const p2 = tableEditor.inputs.find((p) => p.name === 'schedule')!;
    const notRun = customEditorFor(tableEditor, 'schedule')!(p2, context());
    expect(notRun.element.querySelector('button') != null, true, 'not computed yet: a Run up to here button');
    const unconnected = customEditorFor(tableEditor, 'schedule')!(p2, context({isConnected: () => false}));
    expect(unconnected.element.textContent?.includes('Connect table'), true, 'not wired: asks to connect it');
  });

  test('a flow input wired into a list<string> input is a nullable string list that keeps commas in its items',
    async () => {
      const e = makeEditor();
      try {
        const script = DG.Script.create('//name: ffStringListTest\n//language: javascript\n' +
          '//input: list<string> groups {nullable: true}\n//output: int n\nn = (groups ?? []).length;\n');
        const func = await addNode(e.flow, ensureFuncNodeType(script));
        const input = await addNode(e.flow, 'Inputs/Property Input');
        await e.flow.addConnectionByKeys(input.id, 'value', func.id, 'groups');
        expect(input.properties['propertyType'], 'string_list',
          'a list<string> parameter is mimicked as a string list');
        expect(input.properties['nullable'], true, 'and stays nullable');
        input.properties['defaultValue'] = ['Dose-response Curve (ELISA, MSD) · qAC50', 'SPR kinetics · KD'];
        expect(JSON.stringify(resolveInputValue(input).value),
          JSON.stringify(['Dose-response Curve (ELISA, MSD) · qAC50', 'SPR kinetics · KD']),
          'a list value is passed as it is, commas inside items included');
        const line = emitScript(e.flow, SETTINGS, {}).split('\n').find((l) => l.includes(' groups'))!;
        expect(line.startsWith('//input: string_list groups'), true, `a string list in the run dialog (${line})`);
        expect(line.includes('='), false, 'with no default: the platform reads none for a list');
        input.properties['defaultValue'] = [];
        expect(resolveInputValue(input).ok, true, 'an empty list is fine: the input is nullable');
      } finally {
        destroyEditor(e);
      }
    });

  test('a flow input wired into an annotated input is edited with the editor, and keeps the value', async () => {
    if (!editorDeployed()) return;
    const e = makeEditor();
    const panel = new PropertyPanel(e.flow);
    document.body.appendChild(panel.root);
    try {
      const func = await addNode(e.flow, ensureFuncNodeType(annotatedScript(EDITOR)));
      const input = await addNode(e.flow, 'Inputs/Property Input');
      await e.flow.addConnectionByKeys(input.id, 'value', func.id, 'schedule');
      input.properties['defaultValue'] = '0 * * * *';
      panel.showNode(func);
      const row = (): Element | null => panel.root.querySelector('[data-param="schedule"]');
      const select = (): HTMLSelectElement | null => row()?.querySelector('select.pp-cron-preset-select') ?? null;
      expect(await until(() => select() != null, 10_000), true,
        'the connected input shows the editor, not "connected"');
      expect(select()!.value, '0 * * * *', 'with the flow input\'s value');
      select()!.value = '*/15 * * * *';
      select()!.dispatchEvent(new Event('change', {bubbles: true}));
      expect(await until(() => input.properties['defaultValue'] === '*/15 * * * *'), true,
        'the edit is written into the flow input');

      const header = emitScript(e.flow, SETTINGS, {}).split('\n').filter((l) => l.startsWith('//input:'));
      expect(header.some((l) => l.includes('editor')), false, 'the run dialog asks for a plain input');
    } finally {
      panel.root.remove();
      destroyEditor(e);
    }
  });

  test('a flow input wired into an input with an editor shows that editor on its card and in its own panel',
    async () => {
      if (!editorDeployed()) return;
      const container = ui.div([], {style: {width: '1000px', height: '700px', position: 'absolute', left: '-10000px'}});
      document.body.appendChild(container);
      const edits: string[] = [];
      const flow = new FlowEditor(container, {onGraphEdited: (e) => {
        if (e.kind === 'params-changed') edits.push(e.nodeId);
      }});
      const panel = new PropertyPanel(flow);
      document.body.appendChild(panel.root);
      try {
        const func = await addNode(flow, ensureFuncNodeType(annotatedScript(EDITOR)));
        const input = await addNode(flow, 'Inputs/Property Input');
        await flow.addConnectionByKeys(input.id, 'value', func.id, 'schedule');
        input.properties['defaultValue'] = '0 * * * *';
        flow.rebuildValueEditor(input.id);
        const control = input.controls['value'] as InputValueControl;
        const card = control.element()!;
        document.body.appendChild(card);
        expect(card.getAttribute('data-testid'), 'ff-editor-box', 'the card hosts the editor in a resizable box');
        const select = (): HTMLSelectElement | null => card.querySelector('select.pp-cron-preset-select');
        expect(await until(() => select() != null, 10_000), true, 'the function input\'s editor is on the card');
        expect(select()!.value, '0 * * * *', 'showing the flow input\'s value');
        const fits = (): boolean =>
          card.scrollWidth <= card.clientWidth + 1 && card.scrollHeight <= card.clientHeight + 1;
        expect(await until(fits, 3000), true, 'the box holds the whole editor ' +
          `(${card.scrollWidth}×${card.scrollHeight} in ${card.clientWidth}×${card.clientHeight})`);
        expect(input.properties['editorSize'] === undefined, true,
          'a size fitted to the editor is not stored: only a resize is');
        select()!.value = '*/5 * * * *';
        select()!.dispatchEvent(new Event('change', {bubbles: true}));
        expect(await until(() => input.properties['defaultValue'] === '*/5 * * * *'), true,
          'an edit on the card is the flow input\'s value');
        expect(edits.includes(input.id), true, 'reported as an edit of the flow input');

        panel.showNode(input);
        expect(await until(() => panel.root.querySelector('[data-param="Value"] select.pp-cron-preset-select') != null,
          10_000), true, 'the flow input\'s own panel shows the editor too');

        panel.showNode(func);
        const inPanel = (): HTMLSelectElement | null =>
          panel.root.querySelector('[data-param="schedule"] select.pp-cron-preset-select');
        expect(await until(() => inPanel() != null, 10_000), true, 'the function node\'s panel shows it');
        inPanel()!.value = '*/15 * * * *';
        inPanel()!.dispatchEvent(new Event('change', {bubbles: true}));
        expect(await until(() => select()?.value === '*/15 * * * *'), true, 'the card follows an edit in the panel');
        expect(control.element() === card, true, 'brought in step, not rebuilt');

        const wire = flow.getConnections().find((c) => c.source === input.id)!;
        await flow.editor.removeConnection(wire.id);
        const plain = control.element()!;
        expect(plain !== card && plain.querySelector('.pp-cron-preset-select') == null, true,
          'unwired, the card goes back to the plain input');
      } finally {
        panel.root.remove();
        flow.destroy();
        container.remove();
      }
    });

  test('a card editor is rebuilt after a run only when the tables it reads have changed', async () => {
    if (!editorDeployed()) return;
    const container = ui.div([], {style: {width: '1000px', height: '700px', position: 'absolute', left: '-10000px'}});
    document.body.appendChild(container);
    let captured: DG.DataFrame | null = null;
    const flow = new FlowEditor(container, {getUpstreamTable: () => captured});
    try {
      const func = await addNode(flow, ensureFuncNodeType(annotatedScript(EDITOR)));
      const table = await addNode(flow, 'Inputs/Table Input');
      const input = await addNode(flow, 'Inputs/Property Input');
      await flow.addConnectionByKeys(table.id, 'table', func.id, 'table');
      await flow.addConnectionByKeys(input.id, 'value', func.id, 'schedule');
      const control = input.controls['value'] as InputValueControl;
      const first = control.element();
      flow.refreshWiredEditors();
      expect(control.element() === first, true, 'nothing changed: the card keeps its editor');
      captured = DG.DataFrame.fromCsv('a\n1\n');
      expect(control.isStale(), true, 'a captured table it did not have makes it stale');
      flow.refreshWiredEditors();
      const second = control.element();
      expect(second !== first, true, 'rebuilt from the new table');
      flow.refreshWiredEditors();
      expect(control.element() === second, true, 'and kept while that table stays');
    } finally {
      flow.destroy();
      container.remove();
    }
  });
});
