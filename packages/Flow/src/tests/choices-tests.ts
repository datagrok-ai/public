/** Choices of an input node: parsing, checking, the panel field that commits them, the node's editor
 *  rebuilt from them, and the header line the run dialog reads them from. */

import * as DG from 'datagrok-api/dg';
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {literalChoiceItems} from '../utils/choice-refs';
import {registerBuiltinNodes, registerAllFunctions} from '../rete/node-factory';
import {checkChoices, headerChoicesRef} from '../panel/choice-input-processor';
import {PropertyPanel} from '../panel/property-panel';
import {emitScript} from '../compiler/script-emitter';
import {InputValueControl} from '../rete/nodes/input-value-control';
import {FlowNode} from '../rete/scheme';
import {makeEditor, destroyEditor, addNode, until} from './test-utils';

const SETTINGS = {name: 'ChoicesTest', description: '', tags: []};
/** A function returning choices that a stand with Chem has; the tests needing it skip without it. */
const CHOICES_FUNC = 'Chem:getMpoProfileNames';

function choicesFuncDeployed(): boolean {
  return DG.Func.find({package: 'Chem', name: 'getMpoProfileNames'}).length > 0;
}

/** Types into the panel's Choices field and commits it as leaving the field does. */
function typeChoices(panel: PropertyPanel, text: string, how: 'blur' | 'enter' = 'blur'): HTMLInputElement {
  const field = panel.root.querySelector('[data-param="Choices"] input') as HTMLInputElement;
  field.value = text;
  field.dispatchEvent(new Event('input', {bubbles: true}));
  if (how === 'enter')
    field.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', bubbles: true}));
  else
    field.dispatchEvent(new Event('change', {bubbles: true}));
  return field;
}

/** The node card's value editor, built as a render would build it. */
function nodeEditor(node: FlowNode): HTMLElement | null {
  const ctl = node.controls['value'];
  return ctl instanceof InputValueControl ? ctl.element() : null;
}

category('Flow: input choices', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('literal choices: quoted, bracketed or bare items', async () => {
    expect(literalChoiceItems('"AA", "BB"').join('|'), 'AA|BB', 'double-quoted items lose their quotes');
    expect(literalChoiceItems('[\'a\', \'b\']').join('|'), 'a|b', 'a bracketed, single-quoted list');
    expect(literalChoiceItems(' x ,, y ').join('|'), 'x|y', 'blanks dropped, items trimmed');
    expect(literalChoiceItems('').length, 0);
  });

  test('checkChoices: literal lists, function references and what is wrong with them', async () => {
    expect((checkChoices('"AA", "BB"') as {value: string}).value, 'AA, BB', 'a literal list is stored bare');
    expect((checkChoices('') as {value: string}).value, '', 'empty clears the choices');
    expect('error' in checkChoices('NoSuchPackage:noSuchFunction()'), true, 'an unknown function is an error');
    expect((checkChoices('a:b') as {value: string}).value, 'a:b',
      'a colon without parentheses that names no function is a literal value');
    if (!choicesFuncDeployed()) return;
    const typed = checkChoices('chem:getMpoProfileNames()') as {value: string};
    expect(typed.value, `${CHOICES_FUNC}()`, 'spelled with the function\'s real name');
    expect((checkChoices('Chem:getMpoProfileNames') as {value: string}).value, `${CHOICES_FUNC}()`,
      'parentheses are optional');
    expect('error' in checkChoices('Bio:getMpoProfileNames()'), true, 'the function must be in the package named');
    expect('error' in checkChoices('Chem:getMpoProfileNames(a, b)'), true, 'too many arguments');
    expect(headerChoicesRef('hem:getMpoProfileNames('), `${CHOICES_FUNC}()`,
      'a reference adopted from a parameter (first character lost) is spelled whole in the header');
  });

  test('the Choices field commits on leaving it and rebuilds the node\'s editor as a dropdown', async () => {
    const e = makeEditor();
    const panel = new PropertyPanel(e.flow);
    document.body.appendChild(panel.root);
    try {
      const node = await addNode(e.flow, 'Inputs/String Input');
      panel.showNode(node);
      expect(nodeEditor(node)?.querySelector('select') == null, true, 'no choices: a text field on the card');

      const field = panel.root.querySelector('[data-param="Choices"] input') as HTMLInputElement;
      field.value = '"AA", "BB"';
      field.dispatchEvent(new Event('input', {bubbles: true}));
      expect(String(node.properties['choices'] ?? ''), '', 'typing alone stores nothing');

      typeChoices(panel, '"AA", "BB"', 'blur');
      expect(node.properties['choices'], 'AA, BB', 'committed on leaving the field');
      const select = nodeEditor(node)?.querySelector('select');
      expect(!!select, true, 'the card now offers a dropdown');
      // An unset value shows as an empty entry of its own.
      const options = Array.from(select!.querySelectorAll('option')).map((o) => o.textContent).filter(Boolean);
      expect(options.join('|'), 'AA|BB', 'with the committed items, unquoted');
      expect(await until(() => !!panel.root.querySelector('[data-param="Value"] select')), true,
        'the panel\'s Value row is a dropdown too');

      typeChoices(panel, 'x, y, z', 'enter');
      expect(node.properties['choices'], 'x, y, z', 'Enter commits too');
    } finally {
      panel.root.remove();
      destroyEditor(e);
    }
  });

  test('an unknown choices function is flagged on the field and never reaches the node', async () => {
    const e = makeEditor();
    const panel = new PropertyPanel(e.flow);
    document.body.appendChild(panel.root);
    try {
      const node = await addNode(e.flow, 'Inputs/String Input');
      node.properties['choices'] = 'AA, BB';
      panel.showNode(node);
      typeChoices(panel, 'NoSuchPackage:noSuchFunction()');
      expect(node.properties['choices'], 'AA, BB', 'the node keeps its last valid choices');
      const root = panel.root.querySelector('[data-param="Choices"]')!;
      expect(await until(() => root.querySelector('.d4-invalid, .ui-input-invalid, [class*="invalid"]') != null),
        true, 'the field says it is invalid');
    } finally {
      panel.root.remove();
      destroyEditor(e);
    }
  });

  test('a choices function fills the card\'s dropdown and rides the header the run dialog reads', async () => {
    if (!choicesFuncDeployed()) return;
    const e = makeEditor();
    const panel = new PropertyPanel(e.flow);
    document.body.appendChild(panel.root);
    try {
      const node = await addNode(e.flow, 'Inputs/String Input');
      node.properties['paramName'] = 'profile';
      panel.showNode(node);
      typeChoices(panel, 'chem:getMpoProfileNames()');
      expect(node.properties['choices'], `${CHOICES_FUNC}()`, 'stored with the real name');
      const editor = nodeEditor(node)!;
      expect(await until(() => {
        const options = Array.from(editor.querySelectorAll('option')).map((o) => o.textContent ?? '');
        return options.length > 0 && options.every((o) => !o.includes('getMpoProfileNames'));
      }, 20_000), true, 'the card\'s dropdown holds what the function returned');

      const script = emitScript(e.flow, SETTINGS, {});
      expect(script.includes(`//input: string profile {choices: ${CHOICES_FUNC}()}`), true,
        `the header names the function (${script.split('\n').find((l) => l.includes('profile'))})`);
      const parsed = DG.Script.create(script);
      const choices = parsed.inputs.find((p) => p.name === 'profile')?.choices ?? [];
      expect(choices.length === 1 && String(choices[0]).includes('getMpoProfileNames'), true,
        'the platform reads it as a reference, not a literal list');

      node.properties['choices'] = 'NoSuchPackage:noSuchFunction()';
      expect(emitScript(e.flow, SETTINGS, {}).includes('choices'), false, 'an unresolvable reference stays out');
      node.properties['choices'] = '"AA", "BB"';
      expect(emitScript(e.flow, SETTINGS, {}).includes('{choices: ["AA", "BB"]}'), true,
        'a literal list, unquoted once');
    } finally {
      panel.root.remove();
      destroyEditor(e);
    }
  });
});
