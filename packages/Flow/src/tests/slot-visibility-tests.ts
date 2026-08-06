/** Default-hidden primitive input rows and the per-node "Shown inputs"
 *  overrides. The data layer (slots, seeds, pass-throughs, compile) is never
 *  touched by visibility — only which rows the card renders. */
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes, registerAllFunctions, getRegisteredFuncs} from '../rete/node-factory';
import {hiddenSocketRow} from '../rete/scheme';
import {isPrimitiveSlotType} from '../types/type-map';
import {serializeFlow, deserializeFlow} from '../serialization/flow-serializer';
import {tid} from '../utils/test-ids';
import {makeEditor, destroyEditor, addNode, until} from './test-utils';

function funcTypeName(name: string): string | null {
  return getRegisteredFuncs().find((f) => f.func.name === name)?.nodeTypeName ?? null;
}

const SETTINGS = {scriptName: 'test', scriptDescription: '', tags: ['funcflow']};

function inputRow(container: HTMLElement, key: string): boolean {
  return !!container.querySelector(`[data-testid="${tid('socket-input', key)}"]`);
}

function outputRow(container: HTMLElement, key: string): boolean {
  return !!container.querySelector(`[data-testid="${tid('socket-output', key)}"]`);
}

category('Flow: slot visibility', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('isPrimitiveSlotType: scalars yes, structural types no', async () => {
    for (const t of ['string', 'int', 'double', 'num', 'bool', 'datetime'])
      expect(isPrimitiveSlotType(t), true, `${t} is primitive`);
    for (const t of ['dataframe', 'column', 'column_list', 'file', 'list', 'string_list', 'list<string>', 'object'])
      expect(isPrimitiveSlotType(t), false, `${t} is not primitive`);
  });

  test('primitive input rows and their pass-throughs are hidden by default', async () => {
    const typeName = funcTypeName('AddNewColumn');
    if (!typeName) return; // not registered on this server — skip
    const e = makeEditor();
    try {
      const node = await addNode(e.flow, typeName);
      expect(node.inputSlotShownByDefault('table'), true, 'dataframe input shown by default');
      expect(node.inputSlotShownByDefault('name'), false, 'string input hidden by default');
      expect(node.inputSlotShownByDefault('expression'), false, 'string input hidden by default');
      // Data layer untouched — slots, seeds, and pass-throughs all still exist.
      expect('name' in node.inputs && 'name__pt' in node.outputs, true, 'hidden slot stays data-carrying');

      await until(() => inputRow(e.container, 'table'));
      expect(inputRow(e.container, 'name'), false, 'no row for a primitive input');
      expect(inputRow(e.container, 'expression'), false, 'no row for a primitive input');
      expect(outputRow(e.container, 'name__pt'), false, 'no pass-through row either');
      expect(outputRow(e.container, 'table__pt'), true, 'dataframe pass-through still renders');
      expect(outputRow(e.container, 'result'), true, 'real output always renders');
    } finally {
      destroyEditor(e);
    }
  });

  test('setInputShown reveals the row and its pass-through; default resets tidy', async () => {
    const typeName = funcTypeName('AddNewColumn');
    if (!typeName) return;
    const e = makeEditor();
    try {
      const node = await addNode(e.flow, typeName);
      await until(() => inputRow(e.container, 'table'));

      await e.flow.setInputShown(node.id, 'expression', true);
      expect((node.properties['shownSlots'] as Record<string, boolean>)['expression'], true, 'override stored');
      await until(() => inputRow(e.container, 'expression'));
      expect(outputRow(e.container, 'expression__pt'), true, 'pass-through follows the input');
      expect(inputRow(e.container, 'name'), false, 'other primitives stay hidden');

      await e.flow.setInputShown(node.id, 'expression', false);
      expect('shownSlots' in node.properties, false, 'back-to-default override removed entirely');
      await until(() => !inputRow(e.container, 'expression'));
      expect(outputRow(e.container, 'expression__pt'), false, 'pass-through hidden again');

      // Non-primitives can be hidden too — the checkbox works both ways.
      await e.flow.setInputShown(node.id, 'table', false);
      await until(() => !inputRow(e.container, 'table'));
      expect(outputRow(e.container, 'table__pt'), false, 'its pass-through goes with it');
    } finally {
      destroyEditor(e);
    }
  });

  test('a wired primitive input renders regardless of the default', async () => {
    const typeName = funcTypeName('AddNewColumn');
    if (!typeName) return;
    const e = makeEditor();
    try {
      const str = await addNode(e.flow, 'Inputs/String Input', 0, 0);
      const node = await addNode(e.flow, typeName, 300, 0);
      await e.flow.addConnectionByKeys(str.id, 'value', node.id, 'name');
      await until(() => inputRow(e.container, 'name'));
      expect(inputRow(e.container, 'expression'), false, 'unwired primitives stay hidden');
      // The predicate agrees: wired → not hidden, unwired → hidden.
      expect(hiddenSocketRow(node, 'input', 'name', (s, k) => s === 'input' && k === 'name'), false);
      expect(hiddenSocketRow(node, 'input', 'name', () => false), true);
    } finally {
      destroyEditor(e);
    }
  });

  test('shownSlots overrides survive a save/load round-trip', async () => {
    const typeName = funcTypeName('AddNewColumn');
    if (!typeName) return;
    const e = makeEditor();
    try {
      const node = await addNode(e.flow, typeName);
      await e.flow.setInputShown(node.id, 'expression', true);
      const doc = serializeFlow(e.flow, SETTINGS);
      await deserializeFlow(doc, e.flow);
      const re = e.flow.getNodes().find((n) => n.dgTypeName === typeName)!;
      expect(re.inputSlotShown('expression'), true, 'override restored');
      expect(re.inputSlotShown('name'), false, 'default untouched');
      await until(() => inputRow(e.container, 'expression'));
    } finally {
      destroyEditor(e);
    }
  });

  test('duplicating a node carries its shownSlots overrides', async () => {
    const typeName = funcTypeName('AddNewColumn');
    if (!typeName) return;
    const e = makeEditor();
    try {
      const node = await addNode(e.flow, typeName);
      await e.flow.setInputShown(node.id, 'expression', true);
      const copies = await e.flow.duplicateNodes([node.id]);
      expect(copies.length, 1);
      expect(copies[0].inputSlotShown('expression'), true, 'copy shows the same rows');
    } finally {
      destroyEditor(e);
    }
  });

  test('the ⋯ indicator marks hidden toggleable inputs and pops the checkboxes', async () => {
    const typeName = funcTypeName('AddNewColumn');
    if (!typeName) return;
    const e = makeEditor();
    try {
      const node = await addNode(e.flow, typeName);
      const dots = (): HTMLElement | null =>
        e.container.querySelector(`[data-testid="${tid('node-more-inputs')}"]`);
      expect(await until(() => dots() != null), true, 'indicator renders while primitives are hidden');

      dots()!.dispatchEvent(new MouseEvent('click', {bubbles: true, cancelable: true}));
      const checkbox = (): HTMLElement | null =>
        Array.from(document.querySelectorAll<HTMLElement>('.d4-menu-item-label'))
          .find((el) => el.textContent?.trim() === (node.inputs['expression']?.label ?? 'expression')) ?? null;
      expect(await until(() => checkbox() != null), true, 'clicking pops the Shown inputs checkboxes');
      checkbox()!.click();
      expect(await until(() => node.inputSlotShown('expression')), true, 'the checkbox toggles the row');

      // Show every remaining toggleable input — the indicator goes away.
      for (const k of Object.keys(node.inputs)) {
        if (!k.startsWith('__exec') && !node.hiddenInputs.has(k) && !node.inputSlotShown(k))
          await e.flow.setInputShown(node.id, k, true);
      }
      expect(await until(() => dots() == null), true, 'nothing hidden — no indicator');
    } finally {
      for (const el of Array.from(document.querySelectorAll('.d4-menu-popup, .d4-menu-dropdown')))
        el.remove();
      destroyEditor(e);
    }
  });

  test('the ⋯ indicator ignores force-hidden inputs (OpenFile stays clean)', async () => {
    const typeName = funcTypeName('OpenFile');
    if (!typeName) return;
    const e = makeEditor();
    try {
      await addNode(e.flow, typeName);
      await until(() => !!e.container.querySelector('.ff-node'));
      expect(!!e.container.querySelector(`[data-testid="${tid('node-more-inputs')}"]`), false,
        'panel-only inputs are not "hidden inputs" to the user');
    } finally {
      destroyEditor(e);
    }
  });

  test('node context menu offers Shown inputs checkboxes', async () => {
    const typeName = funcTypeName('AddNewColumn');
    if (!typeName) return;
    const e = makeEditor();
    try {
      const node = await addNode(e.flow, typeName);
      await until(() => !!e.container.querySelector('.ff-node'));
      const nodeEl = e.container.querySelector('.ff-node') as HTMLElement;
      nodeEl.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, clientX: 10, clientY: 10}));
      const menuItem = (): HTMLElement | null =>
        Array.from(document.querySelectorAll<HTMLElement>('.d4-menu-item-label'))
          .find((el) => el.textContent?.trim() === 'Shown inputs') ?? null;
      expect(await until(() => menuItem() != null), true, 'Shown inputs group is in the node menu');
      menuItem()!.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      const checkbox = (): HTMLElement | null =>
        Array.from(document.querySelectorAll<HTMLElement>('.d4-menu-item-label'))
          .find((el) => el.textContent?.trim() === (node.inputs['expression']?.label ?? 'expression')) ?? null;
      expect(await until(() => checkbox() != null), true, 'each input gets a checkbox item');
    } finally {
      for (const el of Array.from(document.querySelectorAll('.d4-menu-popup, .d4-menu-dropdown')))
        el.remove();
      destroyEditor(e);
    }
  });
});
