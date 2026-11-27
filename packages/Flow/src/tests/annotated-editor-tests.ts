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
});
