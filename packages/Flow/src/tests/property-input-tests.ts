/** Property Input — the universal input node that MIMICS the function input it
 *  connects to: the target parameter's DG.Property drives the node's type, value
 *  editor, qualifiers, and context panel (`adoptInputProperty`), re-picked on every
 *  connect; non-function targets are configured manually via the panel's Type combo. */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes, registerAllFunctions, createNode, ensureFuncNodeType,
  findNodeTypesProducingOutput} from '../rete/node-factory';
import {PROPERTY_INPUT_TYPE, applyPropertyInputShape, adoptInputProperty} from '../rete/nodes/input-nodes';
import {hostsInlineSketcher, hostsHelmEditor, FlowNode} from '../rete/scheme';
import {TypedSocket} from '../rete/sockets';
import {inputValueProperty, resolveInputValue, buildInputValueEditor} from '../utils/input-values';
import {isChoiceQueryRef} from '../utils/choice-refs';
import {emitScript} from '../compiler/script-emitter';
import {serializeFlow, deserializeFlow} from '../serialization/flow-serializer';
import {PropertyPanel} from '../panel/property-panel';
import {builtinNodeDesc} from '../rete/builtin-catalog';
import {makeEditor, destroyEditor, addNode, until} from './test-utils';

function outSocketType(node: FlowNode): string {
  return (node.outputs['value'] as {socket: TypedSocket} | undefined)?.socket.dgType ?? '?';
}

function funcTypeName(pkg: string, name: string): string {
  const func = DG.Func.find({package: pkg, name})[0];
  if (!func) throw new Error(`${pkg}:${name} is not on this stand`);
  return ensureFuncNodeType(func);
}

category('Flow: property input', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('starts dynamic and leads the reverse suggestions for a function-input drag', async () => {
    const node = createNode(PROPERTY_INPUT_TYPE)!;
    expect(node.dgNodeType, 'input');
    expect(node.dgOutputType, 'dynamic');
    expect(outSocketType(node), 'dynamic');
    expect(String(node.properties['propertyType']), 'dynamic');

    // Dragging out of a DG-function input → Property Input first (it will mimic the param).
    const fromFunc = findNodeTypesProducingOutput('dataframe', {targetIsFuncInput: true});
    expect(fromFunc[0].typeName, PROPERTY_INPUT_TYPE, 'the mimic leads a function-input drag');
    expect(fromFunc[1].typeName, 'Inputs/Table Input', 'the matching input node still follows');

    // A non-function input drag keeps the old order — the exact-matching Input node leads.
    const plain = findNodeTypesProducingOutput('dataframe');
    expect(plain[0].typeName, 'Inputs/Table Input');
    expect(plain.some((c) => c.typeName === PROPERTY_INPUT_TYPE), true,
      'still offered (its dynamic output is compatible), just not first');
  });

  test('adopts the connected function input — the Gasteiger molecule param', async () => {
    const e = makeEditor();
    try {
      const typeName = funcTypeName('Chem', 'ChemistryGasteigerPartialCharges');
      const func = await addNode(e.flow, typeName);
      const input = await addNode(e.flow, PROPERTY_INPUT_TYPE, 400, 0);
      await e.flow.addConnectionByKeys(input.id, 'value', func.id, 'mol');

      expect(input.dgOutputType, 'string', 'adopted the param type');
      expect(outSocketType(input), 'string', 'the output socket followed');
      expect(String(input.properties['semType']), 'Molecule', 'the semType qualifier rode along');
      expect(hostsInlineSketcher(input), true, 'a Molecule string input hosts the inplace sketcher');
      expect(inputValueProperty(input)!.semType, 'Molecule', 'the value editor property is tagged');
      expect(String(input.properties['paramName']), 'mol', 'param name follows the mimicked parameter');
      expect(input.label, 'Mol Input', 'auto-titled after the parameter');
      expect(String(input.properties['defaultValue']).startsWith('COc1'), true,
        'seeded with the declared default molecule');
      expect(resolveInputValue(input).ok, true, 'the seeded value resolves — the node is run-ready');
    } finally {
      destroyEditor(e);
    }
  });

  test('re-picks the property when connected elsewhere; user renames survive', async () => {
    const e = makeEditor();
    try {
      const typeName = funcTypeName('Chem', 'ChemistryGasteigerPartialCharges');
      const func = await addNode(e.flow, typeName);
      const input = await addNode(e.flow, PROPERTY_INPUT_TYPE, 400, 0);
      await e.flow.addConnectionByKeys(input.id, 'value', func.id, 'mol');
      const molWire = e.flow.getConnections().find((c) => c.source === input.id)!;

      // Re-route to the int param → the whole shape re-picks, stale value included.
      await e.flow.editor.removeConnection(molWire.id);
      await e.flow.addConnectionByKeys(input.id, 'value', func.id, 'contours');
      expect(input.dgOutputType, 'int');
      expect(outSocketType(input), 'int');
      expect('semType' in input.properties, false, 'string qualifiers dropped');
      expect('min' in input.properties && 'max' in input.properties, true, 'numeric qualifiers appeared');
      expect(String(input.properties['paramName']), 'contours', 'auto param name followed');
      expect(input.label, 'Contours Input', 'auto title followed');
      expect(String(input.properties['defaultValue']), '10',
        'the SMILES seeded for the previous type was replaced by the int default');

      // A user rename sticks across the next re-adoption.
      input.properties['paramName'] = 'myCount';
      input.label = 'My Special Input';
      const wire = e.flow.getConnections().find((c) => c.source === input.id)!;
      await e.flow.editor.removeConnection(wire.id);
      await e.flow.addConnectionByKeys(input.id, 'value', func.id, 'mol');
      expect(input.dgOutputType, 'string', 'the shape still re-picked');
      expect(String(input.properties['paramName']), 'myCount', 'the renamed param name survived');
      expect(input.label, 'My Special Input', 'the renamed title survived');
    } finally {
      destroyEditor(e);
    }
  });

  test('a Macromolecule param routes the node to the Helm editor', async () => {
    const e = makeEditor();
    try {
      const typeName = funcTypeName('Bio', 'toAtomicLevelSingleSeq');
      const func = await addNode(e.flow, typeName);
      const input = await addNode(e.flow, PROPERTY_INPUT_TYPE, 400, 0);
      await e.flow.addConnectionByKeys(input.id, 'value', func.id, 'sequence');
      expect(input.dgOutputType, 'string');
      expect(String(input.properties['semType']), 'Macromolecule');
      expect(hostsHelmEditor(input), true, 'the node body hosts the HELM editor box');
    } finally {
      destroyEditor(e);
    }
  });

  test('a non-function target leaves the node for manual panel configuration', async () => {
    const e = makeEditor();
    try {
      const input = await addNode(e.flow, PROPERTY_INPUT_TYPE);
      const toStr = await addNode(e.flow, 'Utilities/ToString', 400, 0);
      await e.flow.addConnectionByKeys(input.id, 'value', toStr.id, 'value');
      expect(input.dgOutputType, 'dynamic', 'nothing to mimic on a utility');
      expect(String(input.properties['paramName']), 'value', 'untouched');

      // What the panel's Type combo does: set the type, re-derive the shape.
      input.properties['propertyType'] = 'double';
      applyPropertyInputShape(input);
      expect(input.dgOutputType, 'double');
      expect(outSocketType(input), 'double');
      expect('min' in input.properties && 'showSlider' in input.properties, true,
        'numeric qualifier rows appear in the panel');
      expect('choices' in input.properties, false, 'string qualifiers are not offered');
    } finally {
      destroyEditor(e);
    }
  });

  test('the context panel offers Type plus exactly the adopted family\'s qualifiers', async () => {
    const e = makeEditor();
    try {
      const typeName = funcTypeName('Chem', 'ChemistryGasteigerPartialCharges');
      const func = await addNode(e.flow, typeName);
      const molInput = await addNode(e.flow, PROPERTY_INPUT_TYPE, 400, 0);
      await e.flow.addConnectionByKeys(molInput.id, 'value', func.id, 'mol');
      const intInput = await addNode(e.flow, PROPERTY_INPUT_TYPE, 400, 200);
      await e.flow.addConnectionByKeys(intInput.id, 'value', func.id, 'contours');

      const panel = new PropertyPanel(e.flow);
      const row = (label: string): HTMLElement | null =>
        panel.root.querySelector(`[data-param="${label}"]`);

      panel.showNode(molInput);
      expect(row('Type') != null, true, 'the mimicked-type combo renders');
      expect(row('SemType') != null, true, 'string family: SemType row');
      expect(row('Choices (comma-sep)') != null, true, 'string family: Choices row');
      expect(row('Min'), null, 'no numeric rows on a string');

      panel.showNode(intInput);
      expect(row('Min') != null && row('Max') != null, true, 'int family: Min/Max rows');
      expect(row('Show Slider') != null, true, 'int family: slider toggle');
      expect(row('SemType'), null, 'no string rows on an int');

      // The combo is Property Input-specific — a plain String Input has no Type row.
      const plain = await addNode(e.flow, 'Inputs/String Input', 0, 400);
      panel.showNode(plain);
      expect(row('Type'), null, 'dedicated input nodes keep their fixed type');
    } finally {
      destroyEditor(e);
    }
  });

  test('compiles to the mimicked header lines and the platform parses them', async () => {
    const e = makeEditor();
    try {
      const typeName = funcTypeName('Chem', 'ChemistryGasteigerPartialCharges');
      const func = await addNode(e.flow, typeName);
      const molInput = await addNode(e.flow, PROPERTY_INPUT_TYPE, 400, 0);
      await e.flow.addConnectionByKeys(molInput.id, 'value', func.id, 'mol');
      const intInput = await addNode(e.flow, PROPERTY_INPUT_TYPE, 400, 200);
      await e.flow.addConnectionByKeys(intInput.id, 'value', func.id, 'contours');

      const script = emitScript(e.flow, {name: 'PropInputE2E', description: 'test', tags: ['funcflow']});
      const molLine = script.split('\n').find((l) => l.startsWith('//input: string mol'));
      expect(molLine != null, true, `molecule input line emitted; got:\n${script.split('\n', 8).join('\n')}`);
      expect(molLine!.includes('semType: Molecule'), true, 'the adopted qualifier is in the header');
      expect(script.split('\n').some((l) => l.startsWith('//input: int contours = 10')), true,
        'the int input line carries the adopted default');
      expect(script.includes('ChemistryGasteigerPartialCharges'), true, 'the call is in the body');

      const parsed = DG.Script.create(script);
      const mol = parsed.inputs.find((p) => p.name === 'mol');
      expect(mol != null, true, 'the platform parsed the mimicked parameter');
      expect(mol!.semType, 'Molecule', 'semType survived DG.Script.create');
      expect(String(parsed.inputs.find((p) => p.name === 'contours')?.propertyType), 'int');
    } finally {
      destroyEditor(e);
    }
  }, {timeout: 30000});

  test('a query input\'s choices reference is adopted and resolves to the live list', async () => {
    // Biologics' Assays by Organism: `organism` declares
    // `choices: query("SELECT distinct name FROM biologics.target_organisms")` —
    // a REFERENCE that must ride along (with the query's connection) and resolve
    // in the value editor, exactly like the func node's own panel does.
    const q = await grok.dapi.queries.filter('name = "assaysByOrganism"').first();
    if (q == null) {
      console.warn('Flow: property input: Biologics assaysByOrganism is not on this stand — ' +
        'query-choices test skipped');
      expect(true, true, 'skipped: query unavailable');
      return;
    }
    const e = makeEditor();
    try {
      const func = await addNode(e.flow, ensureFuncNodeType(q));
      const input = await addNode(e.flow, PROPERTY_INPUT_TYPE, 400, 0);
      await e.flow.addConnectionByKeys(input.id, 'value', func.id, 'organism');

      const storedChoices = String(input.properties['choices'] ?? '');
      expect(isChoiceQueryRef(storedChoices), true,
        `the reference is stored verbatim, not mangled into a literal list; got "${storedChoices}"`);
      expect(String(input.properties['choicesConnection'] ?? '') !== '', true,
        'the owning query\'s connection is captured for the resolver');

      let changed = 0;
      const ed = buildInputValueEditor(input, () => changed++)!;
      expect(ed.input instanceof DG.ChoiceInput, true, 'a choices-carrying string edits as a combo');
      const items = (): string[] => {
        try {
          return ((ed.input as DG.ChoiceInput<string>).items ?? []).map((x) => String(x));
        } catch {
          return [];
        }
      };
      const resolved = await until(() =>
        items().length > 0 && !items()[0].toLowerCase().includes('select'), 20000);
      expect(resolved, true, `the reference resolved into real organisms; items: ${items().join(', ')}`);
      expect(changed, 0, 'resolving the list is not a user edit');

      const script = emitScript(e.flow, {name: 'ChoicesRef', description: 'test', tags: ['funcflow']});
      const line = script.split('\n').find((l) => l.startsWith('//input: string organism'))!;
      expect(line != null, true, 'the mimicked input line emitted');
      expect(line.includes('choices:'), false, 'a reference never rides the header');
      expect(DG.Script.create(script).inputs.some((p) => p.name === 'organism'), true,
        'the platform parses the emitted header');
    } finally {
      destroyEditor(e);
    }
  }, {timeout: 60000});

  test('a JSON-escaped query reference (as the platform hands it) still adopts and resolves', async () => {
    // On a live stand the choices can arrive JSON-escaped — `uery(\"SELECT …\"` —
    // and an escape-blind detector reads that as a literal list: the reference is
    // stored mangled and the connection is never captured, so the combo stays bogus.
    const q = await grok.dapi.queries.filter('name = "assaysByOrganism"').first();
    if (q == null) {
      console.warn('Flow: property input: Biologics assaysByOrganism is not on this stand — ' +
        'escaped-choices test skipped');
      expect(true, true, 'skipped: query unavailable');
      return;
    }
    const e = makeEditor();
    try {
      const input = await addNode(e.flow, PROPERTY_INPUT_TYPE);
      const prop = DG.Property.fromOptions({name: 'organism', type: DG.TYPE.STRING} as never);
      prop.choices = ['uery(\\"SELECT distinct name FROM biologics.target_organisms\\"'];
      adoptInputProperty(input, prop, q);

      const stored = String(input.properties['choices'] ?? '');
      expect(stored.includes('\\"'), false, `stored unescaped; got "${stored}"`);
      expect(isChoiceQueryRef(stored), true, 'detected as a query reference');
      expect(String(input.properties['choicesConnection'] ?? '') !== '', true,
        'the owning query\'s connection is captured despite the escaping');

      let changed = 0;
      const ed = buildInputValueEditor(input, () => changed++)!;
      const items = (): string[] => {
        try {
          return ((ed.input as DG.ChoiceInput<string>).items ?? []).map((x) => String(x));
        } catch {
          return [];
        }
      };
      const resolved = await until(() =>
        items().length > 0 && !items()[0].toLowerCase().includes('select'), 20000);
      expect(resolved, true, `the escaped reference resolved into real organisms; items: ${items().join(', ')}`);
      expect(changed, 0, 'resolving the list is not a user edit');
    } finally {
      destroyEditor(e);
    }
  }, {timeout: 60000});

  test('a func-call choices reference on ANY input node resolves in the value editor', async () => {
    // The core mangles a single-entry choices ref (first char + trailing `)` stripped) —
    // the editor must resolve `hem:getMpoProfileNames(` through the registered function.
    const node = createNode('Inputs/String Input')!;
    node.properties['choices'] = 'hem:getMpoProfileNames(';
    let changed = 0;
    const ed = buildInputValueEditor(node, () => changed++)!;
    expect(ed.input instanceof DG.ChoiceInput, true, 'edits as a combo');
    const items = (): string[] => {
      try {
        return ((ed.input as DG.ChoiceInput<string>).items ?? []).map((x) => String(x));
      } catch {
        return [];
      }
    };
    const resolved = await until(() =>
      items().length > 0 && !items()[0].includes('getMpoProfileNames'), 20000);
    expect(resolved, true, `the reference resolved to the function's result; items: ${items().join(', ')}`);
    expect(changed, 0, 'resolving the list is not a user edit');
  }, {timeout: 30000});

  test('save/load and duplicate keep the adopted shape', async () => {
    const e = makeEditor();
    try {
      const typeName = funcTypeName('Chem', 'ChemistryGasteigerPartialCharges');
      const func = await addNode(e.flow, typeName);
      const input = await addNode(e.flow, PROPERTY_INPUT_TYPE, 400, 0);
      await e.flow.addConnectionByKeys(input.id, 'value', func.id, 'mol');

      const doc = serializeFlow(e.flow, {scriptName: 'PropRoundtrip', scriptDescription: '', tags: []});
      const e2 = makeEditor();
      try {
        await deserializeFlow(doc, e2.flow);
        const loaded = e2.flow.getNodes().find((n) => n.dgTypeName === PROPERTY_INPUT_TYPE)!;
        expect(loaded.dgOutputType, 'string', 'type restored from properties');
        expect(outSocketType(loaded), 'string', 'socket restored');
        expect(hostsInlineSketcher(loaded), true, 'semType routing restored');
        expect(e2.flow.getConnections().length, 1, 'the wire survived');
      } finally {
        destroyEditor(e2);
      }

      const [copy] = await e.flow.duplicateNodes([input.id]);
      expect(copy.dgOutputType, 'string', 'a duplicate keeps the adopted type');
      expect(outSocketType(copy), 'string');
    } finally {
      destroyEditor(e);
    }
  }, {timeout: 30000});
});

category('Flow: to semantic value', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('node shape: dynamic in, semantic_value out, Molecule by default', async () => {
    const node = createNode('Utilities/To Semantic Value')!;
    expect(node.dgNodeType, 'utility');
    expect((node.inputs['value'] as {socket: TypedSocket}).socket.dgType, 'dynamic');
    expect((node.outputs['semanticValue'] as {socket: TypedSocket}).socket.dgType, 'semantic_value');
    expect(String(node.properties['semType']), 'Molecule');
    expect(node.requiredProps.includes('semType'), true, 'a blank semType reads as unready');
    expect(builtinNodeDesc('Utilities/To Semantic Value').length > 0, true, 'catalog entry present');
    // The API contract the emission relies on.
    const sv = DG.SemanticValue.fromValueType('CCO', 'Molecule');
    expect(sv.semType, 'Molecule');
    expect(String(sv.value), 'CCO');
  });

  test('emits DG.SemanticValue.fromValueType with the configured semType and runs', async () => {
    const e = makeEditor();
    try {
      const c = await addNode(e.flow, 'Constants/String');
      c.properties['value'] = 'CCO';
      const conv = await addNode(e.flow, 'Utilities/To Semantic Value', 300, 0);
      await e.flow.addConnectionByKeys(c.id, 'value', conv.id, 'value');
      const out = await addNode(e.flow, 'Outputs/Value Output', 600, 0);
      await e.flow.addConnectionByKeys(conv.id, 'semanticValue', out.id, 'value');
      expect(String(out.properties['outputType']), 'semantic_value', 'the Value Output auto-typed');

      const script = emitScript(e.flow, {name: 'SemValE2E', description: 'test', tags: ['funcflow']});
      const convLine = script.split('\n').find((l) => l.includes('DG.SemanticValue.fromValueType('));
      expect(convLine != null, true, `the conversion is in the body:\n${script}`);
      expect(convLine!.includes('"Molecule"'), true, 'the configured semType is passed');

      conv.properties['semType'] = 'Macromolecule';
      const script2 = emitScript(e.flow, {name: 'SemValE2E', description: 'test', tags: ['funcflow']});
      expect(script2.includes('"Macromolecule"'), true, 'the semType property drives the emission');

      const fc = DG.Script.create(script).prepare({});
      await fc.call(undefined, undefined, {processed: true});
      const result: unknown = fc.outputs['result'];
      expect(result != null, true, 'the semantic value flowed to the output');
      if (result instanceof DG.SemanticValue) {
        expect(result.semType, 'Molecule');
        expect(String(result.value), 'CCO');
      }
    } finally {
      destroyEditor(e);
    }
  }, {timeout: 30000});
});
