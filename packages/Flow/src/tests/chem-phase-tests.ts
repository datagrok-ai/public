/** Node contracts the Chem catalog entries must honour on a canvas: sockets match
 *  the panel form, columns resolve against a table, choices references resolve. */

import * as DG from 'datagrok-api/dg';
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes, registerAllFunctions, getRegisteredFuncs, createNode} from '../rete/node-factory';
import {findNodeTypesProducingOutput} from '../rete/node-factory';
import {PropertyPanel} from '../panel/property-panel';
import {
  parseChoiceFuncRef, resolveChoiceFunc, choiceValuesFrom,
} from '../panel/choice-input-processor';
import {isLiteralChoiceList} from '../utils/choice-refs';
import {impliedChoiceDefault} from '../rete/nodes/func-node';
import {
  guessColumnFor, autoMap, parseMapping, serializeMapping, pruneMapping, mappableColumns,
  unmappedProperties, cacheProfileProperties,
} from '../panel/editors/mpo-mapping-editor';
import {
  HIDDEN_FUNC_INPUTS, HIDDEN_FUNC_OUTPUTS, FUNC_WRAPPERS, CUSTOM_FUNC_INPUT_EDITORS,
  hiddenOutputsOf, effectiveFuncInputs, funcWrapperOf,
} from '../utils/func-input-overrides';
import {INCLUDED_FUNC_NQNAMES} from '../rete/included-funcs';
import {inputValueProperty, buildInputValueEditor, resolveInputValue} from '../utils/input-values';
import {missingRequiredInputs, nodeMissingRequirements, FlowNode,
  hostsInlineSketcher} from '../rete/scheme';
import {estimateNodeWidth} from '../rete/graph-layout';
import type {FlowEditor} from '../rete/flow-editor';
import {tid} from '../utils/test-ids';
import {makeEditor, destroyEditor, addNode, until} from './test-utils';

/** Null when the package isn't deployed on this server (tests that need it skip). */
function typeNameOf(nqName: string): string | null {
  return getRegisteredFuncs().find((f) => {
    try {
      return f.func.nqName === nqName;
    } catch {
      return false;
    }
  })?.nodeTypeName ?? null;
}

function inputNames(func: DG.Func): string[] {
  return effectiveFuncInputs(func).map((p) => p.name);
}

/** `missingRequiredInputs` reports labels, not keys — asserting on labels would
 *  silently pass for any caption-carrying parameter. */
function missingKeys(node: FlowNode, flow: FlowEditor): string[] {
  const labels = new Set(missingRequiredInputs(node, (k) => flow.isInputConnected(node.id, k)));
  return node.requiredInputs.filter((k) =>
    labels.has((node.inputs as Record<string, {label?: string} | undefined>)[k]?.label ?? k));
}

category('Flow: choices from a function', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  // Core strips the first character and trailing `)` of a one-entry `choices`
  // string — the damaged-looking fixtures below are deliberate.
  test('parseChoiceFuncRef reads a truncated function reference', async () => {
    const ref = parseChoiceFuncRef('hem:getMpoProfileNames(');
    expect(ref !== null, true, 'parsed');
    expect(ref!.funcName, 'getMpoProfileNames');
    expect(ref!.packagePrefix, 'hem', 'leading character is gone — package is only a hint');
    expect(ref!.args.length, 0);
  });

  test('parseChoiceFuncRef reads literal arguments, with or without the closing paren', async () => {
    const withParen = parseChoiceFuncRef("dmetica:getModels('Absorption')");
    expect(withParen!.funcName, 'getModels');
    expect(withParen!.args.join(','), 'Absorption');
    const truncated = parseChoiceFuncRef("dmetica:getModels('Absorption'");
    expect(truncated!.args.join(','), 'Absorption', 'closing paren is optional');
    const two = parseChoiceFuncRef('kg:f("a", "b"');
    expect(two!.args.join('|'), 'a|b', 'comma-separated, unquoted');
  });

  test('parseChoiceFuncRef rejects everything that is not a function reference', async () => {
    for (const bad of ['', 'Morgan', 'uery("select 1"', 'no-colon(', ':leadingColon(', 'pkg:(', 'pkg:not a name('])
      expect(parseChoiceFuncRef(bad), null, `"${bad}" is not a function reference`);
  });

  test('choiceValuesFrom accepts a list, a dataframe, or a column', async () => {
    expect(choiceValuesFrom(['a', 'b']).join(','), 'a,b');
    expect(choiceValuesFrom([]).length, 0);
    expect(choiceValuesFrom(null).length, 0, 'a null result is not a crash');
    const col = DG.Column.fromStrings('c', ['x', 'y', 'x']);
    expect(choiceValuesFrom(DG.DataFrame.fromColumns([col])).sort().join(','), 'x,y', 'categories of the first column');
    expect(choiceValuesFrom(col).sort().join(','), 'x,y');
  });

  test('resolveChoiceFunc finds the function by name (needs a live catalog)', async () => {
    const ref = parseChoiceFuncRef('hem:getMpoProfileNames(')!;
    const f = resolveChoiceFunc(ref);
    if (!f) return; // Chem not deployed on this stand
    expect(f.name, 'getMpoProfileNames');
  });

  test('a choices-by-function input ends up with real items', async () => {
    const typeName = typeNameOf('Chem:mpoScoreByProfile');
    if (!typeName) return;
    const e = makeEditor();
    const panel = new PropertyPanel(e.flow);
    document.body.appendChild(panel.root);
    try {
      const node = await addNode(e.flow, typeName);
      panel.showNode(node);
      const row = panel.root.querySelector('[data-param="profileName"]');
      expect(!!row, true, 'the Profile row renders');
      // resolves async against server files — slow on a cold cache
      const resolved = await until(() => {
        const opts = Array.from(row!.querySelectorAll('option')).map((o) => o.textContent ?? '');
        return opts.every((o) => !o.includes('getMpoProfileNames'));
      }, 20_000);
      expect(resolved, true, 'the reference was replaced by the function\'s result');
      const options = Array.from(row!.querySelectorAll('option')).map((o) => o.textContent ?? '');
      for (const o of options)
        expect(o.includes('('), false, `"${o}" is a resolved profile name, not the reference itself`);
    } finally {
      panel.root.remove();
      destroyEditor(e);
    }
  });
});

category('Flow: choice defaults', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('impliedChoiceDefault takes the first choice of a required choice input', async () => {
    const required = DG.Property.fromOptions({
      name: 'aggregation', type: 'string', nullable: false, choices: ['Average', 'Sum'],
    } as never);
    expect(impliedChoiceDefault(required), 'Average');
  });

  // `expect`'s expected arg defaults to `true`, so an absent default is compared explicitly.
  test('impliedChoiceDefault leaves optional, choice-less and non-string inputs alone', async () => {
    const nullable = DG.Property.fromOptions({
      name: 'a', type: 'string', nullable: true, choices: ['Average', 'Sum'],
    } as never);
    expect(impliedChoiceDefault(nullable) === undefined, true, 'a nullable choice shows a blank option first');
    const noChoices = DG.Property.fromOptions({name: 'b', type: 'string', nullable: false} as never);
    expect(impliedChoiceDefault(noChoices) === undefined, true, 'a free-text field has nothing to imply');
    const notString = DG.Property.fromOptions(
      {name: 'c', type: 'int', nullable: false, choices: ['1', '2']} as never);
    expect(impliedChoiceDefault(notString) === undefined, true, 'only string inputs render a choice combo');
  });

  test('impliedChoiceDefault never seeds an unresolved choices reference', async () => {
    for (const ref of ['hem:getMpoProfileNames(', 'uery("select name from t"']) {
      const prop = DG.Property.fromOptions(
        {name: 'p', type: 'string', nullable: false, choices: [ref]} as never);
      expect(impliedChoiceDefault(prop) === undefined, true, `${ref} is a reference, not a value`);
    }
    expect(isLiteralChoiceList(['Average', 'Sum']), true);
    expect(isLiteralChoiceList(['hem:getMpoProfileNames(']), false);
    expect(isLiteralChoiceList([]), false);
  });

  test('no fresh node reports a choice parameter as missing', async () => {
    const withChoices = getRegisteredFuncs().filter((f) => {
      try {
        return effectiveFuncInputs(f.func).some((p) =>
          String(p.propertyType) === 'string' && impliedChoiceDefault(p) !== undefined);
      } catch {
        return false;
      }
    }).slice(0, 12); // a sample is enough; building nodes is not free
    expect(withChoices.length > 0, true, 'the catalog has required choice inputs to check');
    for (const info of withChoices) {
      const e = makeEditor();
      try {
        const node = await addNode(e.flow, info.nodeTypeName);
        const missing = missingKeys(node, e.flow);
        for (const p of effectiveFuncInputs(info.func)) {
          if (impliedChoiceDefault(p) === undefined) continue;
          expect(String(node.inputValues[p.name] ?? '').length > 0, true,
            `${info.nodeTypeName}: ${p.name} is seeded with what its combo shows`);
          expect(missing.includes(p.name), false, `${info.nodeTypeName}: ${p.name} is not reported missing`);
        }
      } finally {
        destroyEditor(e);
      }
    }
  });

  test('MPO Score by Profile: aggregation is set, the profile is genuinely unset', async () => {
    const typeName = typeNameOf('Chem:mpoScoreByProfile');
    if (!typeName) return;
    const e = makeEditor();
    try {
      const node = await addNode(e.flow, typeName);
      expect(String(node.inputValues['aggregation'] ?? ''), 'Average', 'seeded with what the combo shows');
      const missing = missingKeys(node, e.flow);
      expect(missing.includes('aggregation'), false, 'so it is not reported missing on a fresh node');
      expect(missing.includes('profileName'), true, 'an unresolved-choices input stays genuinely unset');
    } finally {
      destroyEditor(e);
    }
  });
});

category('Flow: mpo column mapping', () => {
  test('guessColumnFor: exact, then normalized, then containment', async () => {
    expect(guessColumnFor('MW', ['MW', 'MolWeight']), 'MW', 'exact wins');
    expect(guessColumnFor('Molecular Weight', ['molecular_weight']), 'molecular_weight', 'punctuation and case ignored');
    expect(guessColumnFor('LogP', ['Calculated LogP', 'x']), 'Calculated LogP', 'containment');
    expect(guessColumnFor('LogP', ['LogP error', 'LogP']), 'LogP', 'an exact match beats a longer containing one');
    expect(guessColumnFor('LogP', ['Calculated LogP value', 'LogP est']), 'LogP est', 'shortest containing name wins');
    expect(guessColumnFor('MW', ['Activity', 'IC50']), null, 'nothing close enough');
    expect(guessColumnFor('', ['MW']), null);
  });

  test('autoMap fills only blanks, never reusing a column', async () => {
    const props = ['MW', 'LogP', 'TPSA'];
    const cols = ['MW', 'Calculated LogP', 'Activity'];
    const fresh = autoMap(props, cols, {});
    expect(fresh['MW'], 'MW');
    expect(fresh['LogP'], 'Calculated LogP');
    expect('TPSA' in fresh, false, 'nothing matched, so it stays unmapped');

    const decided = autoMap(props, cols, {MW: 'Activity', LogP: ''});
    expect(decided['MW'], 'Activity', 'user mapping kept');
    expect(decided['LogP'], '', 'explicit "unmapped" kept');
  });

  test('mapping round-trips through the stored JSON, and blanks are dropped', async () => {
    expect(serializeMapping({}), '');
    expect(serializeMapping({MW: '', LogP: ''}), '', 'all blank stores nothing');
    const json = serializeMapping({MW: 'Molecular Weight', LogP: ''});
    expect(json, '{"MW":"Molecular Weight"}', 'unmapped properties are omitted');
    expect(parseMapping(json)['MW'], 'Molecular Weight');
    for (const bad of ['', '   ', '{not json', '[1,2]', null, undefined])
      expect(Object.keys(parseMapping(bad)).length, 0, `"${bad}" parses to an empty mapping`);
    expect(parseMapping('{"MW":null}')['MW'], '', 'a null mapping reads as unmapped');
  });

  test('pruneMapping drops entries belonging to a previously chosen profile', async () => {
    const pruned = pruneMapping({MW: 'Molecular Weight', OldProp: 'x'}, ['MW', 'LogP']);
    expect(Object.keys(pruned).join(','), 'MW', 'only current properties survive');
    expect(Object.keys(pruneMapping({}, ['MW'])).length, 0);
  });

  test('only numeric columns are offered', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromFloat32Array('MW', new Float32Array([1, 2])),
      DG.Column.fromStrings('Name', ['a', 'b']),
      DG.Column.fromInt32Array('Count', new Int32Array([1, 2])),
    ]);
    expect(mappableColumns(Array.from(df.columns)).sort().join(','), 'Count,MW');
    expect(mappableColumns(null).length, 0, 'no table → nothing mappable');
    expect(mappableColumns([]).length, 0);
  });

  /** Stand-in for the panel context, so the editor's states are testable without a live node. */
  function fakeCtx(profile: string, columns: DG.Column[] | null) {
    const watchers: Array<(v: unknown) => void> = [];
    const state = {profile};
    return {
      ctx: {
        inputValue: (name: string) => (name === 'profileName' ? state.profile : undefined),
        columns: () => columns,
        watch: (_name: string, cb: (v: unknown) => void) => {watchers.push(cb);},
      },
      setProfile: (p: string) => {
        state.profile = p;
        for (const cb of watchers) cb(p);
      },
    };
  }

  const mappingParam = (): DG.Property =>
    DG.Property.fromOptions({name: 'columnMapping', type: 'string'} as never);

  test('with no profile or no table the editor blocks instead of showing an unfillable form', async () => {
    const factory = CUSTOM_FUNC_INPUT_EDITORS['Chem:mpoScoreByProfile']?.columnMapping;
    expect(!!factory, true, 'registered for the columnMapping input');

    const noProfile = fakeCtx('', [DG.Column.fromFloat32Array('MW', new Float32Array([1]))]);
    const a = factory(mappingParam(), noProfile.ctx);
    a.setValue('');
    await until(() => a.element.getAttribute('data-blocked') === 'true');
    expect(a.element.getAttribute('data-blocked'), 'true', 'the section is blocked');
    expect(a.element.querySelectorAll('[data-mpo-property]').length, 0, 'and renders no empty combos');
    expect((a.element.textContent ?? '').toLowerCase().includes('profile'), true, 'the reason is on screen');

    // a nonexistent profile scores nothing — its own blocked state, reported before the table
    const noSuchProfile = fakeCtx('no such profile', null);
    const b = factory(mappingParam(), noSuchProfile.ctx);
    b.setValue('');
    await until(() => b.element.getAttribute('data-blocked') === 'true');
    expect((b.element.textContent ?? '').includes('no properties'), true, 'says the profile scores nothing');

    const names = await (async (): Promise<string[]> => {
      const f = DG.Func.find({package: 'Chem', name: 'getMpoProfileNames'})[0];
      return f ? ((await f.apply({})) as string[]) : [];
    })();
    if (names.length === 0) return; // no saved profiles on this stand
    const noTable = fakeCtx(names[0], null);
    const c = factory(mappingParam(), noTable.ctx);
    c.setValue('');
    await until(() => (c.element.textContent ?? '').toLowerCase().includes('table'), 20_000);
    expect(c.element.getAttribute('data-blocked'), 'true', 'the section is blocked');
    expect(c.element.querySelectorAll('[data-mpo-property]').length, 0, 'no combos without columns');
    expect((c.element.textContent ?? '').toLowerCase().includes('table'), true, 'says a table is needed');
  });

  test('unmappedProperties reports every property left without a column', async () => {
    const props = ['MW', 'LogP', 'TPSA'];
    expect(unmappedProperties(props, '{"MW":"a","LogP":"b","TPSA":"c"}').length, 0, 'complete');
    expect(unmappedProperties(props, '{"MW":"a"}').join(','), 'LogP,TPSA');
    expect(unmappedProperties(props, '{"MW":"a","LogP":"  "}').join(','), 'LogP,TPSA', 'blank is unmapped');
    expect(unmappedProperties(props, '').join(','), 'MW,LogP,TPSA', 'nothing stored → nothing mapped');
    expect(unmappedProperties([], '').length, 0, 'no properties, nothing to require');
  });

  test('an incomplete mapping makes the node unready, a complete one clears it', async () => {
    const typeName = typeNameOf('Chem:mpoScoreByProfile');
    if (!typeName) return;
    const e = makeEditor();
    try {
      const node = await addNode(e.flow, typeName);
      const isConnected = (k: string): boolean => e.flow.isInputConnected(node.id, k);

      node.inputValues['profileName'] = 'Some profile';
      // Assert only that a mapping requirement is reported, not how many — a blank
      // `columnMapping` can also surface under its own label.
      const unresolved = nodeMissingRequirements(node, isConnected).filter((m) => m.includes('Column mapping'));
      expect(unresolved.length > 0, true, 'an unverifiable mapping blocks the node');

      cacheProfileProperties(node, 'Some profile', ['MW', 'LogP']);
      node.inputValues['columnMapping'] = '{"MW":"Molecular Weight"}';
      const partial = nodeMissingRequirements(node, isConnected).filter((m) => m.includes('Column mapping'));
      expect(partial.length > 0, true, 'a partly-mapped profile is a missing requirement');
      expect(partial.some((m) => m.includes('LogP')), true, 'and it names what is unmapped');

      node.inputValues['columnMapping'] = '{"MW":"Molecular Weight","LogP":"cLogP"}';
      expect(nodeMissingRequirements(node, isConnected).some((m) => m.includes('Column mapping')), false,
        'a complete mapping clears it');

      node.inputValues['profileName'] = 'Another profile';
      expect(nodeMissingRequirements(node, isConnected).some((m) => m.includes('Column mapping')), true,
        'a stale list for a different profile does not count as known');

      // with no profile the profile input is already reported missing — no second complaint
      node.inputValues['profileName'] = '';
      expect(nodeMissingRequirements(node, isConnected).some((m) => m.includes('Column mapping')), false);
    } finally {
      destroyEditor(e);
    }
  });

  test('an unconfigured MPO node is excluded from the run set', async () => {
    const typeName = typeNameOf('Chem:mpoScoreByProfile');
    if (!typeName) return;
    const names = await (async (): Promise<string[]> => {
      const f = DG.Func.find({package: 'Chem', name: 'getMpoProfileNames'})[0];
      return f ? ((await f.apply({})) as string[]) : [];
    })();
    if (names.length === 0) return; // no saved profiles on this stand — skip
    const e = makeEditor();
    try {
      const node = await addNode(e.flow, typeName);
      const isConnected = (k: string): boolean => e.flow.isInputConnected(node.id, k);
      node.inputValues['profileName'] = names[0];
      node.inputValues['columnMapping'] = '';

      expect(nodeMissingRequirements(node, isConnected).some((m) => m.includes('Column mapping')), true,
        'blocked while the profile is still being read');

      // Match across ALL mapping messages — a blank `columnMapping` also surfaces under
      // its own label, and a profile that scores nothing needs no mapping.
      const mappingHints = (): string[] =>
        nodeMissingRequirements(node, isConnected).filter((m) => m.includes('Column mapping'));
      await until(() => mappingHints().some((m) => m.includes('unmapped:')), 20_000);
      const properties = await DG.Func.find({package: 'Chem', name: 'getMpoProfileProperties'})[0]
        .apply({profileName: names[0]}) as string[];
      if (properties.length > 0) {
        const hints = mappingHints();
        expect(hints.length > 0, true, 'an empty mapping never satisfies a profile that scores properties');
        expect(hints.some((m) => m.includes(properties[0])), true, 'and the hint names what is unmapped');
      }
    } finally {
      destroyEditor(e);
    }
  });

  test('changing the profile rebuilds the mapping rows without re-selecting the node', async () => {
    const factory = CUSTOM_FUNC_INPUT_EDITORS['Chem:mpoScoreByProfile']?.columnMapping;
    const names = await (async (): Promise<string[]> => {
      const f = DG.Func.find({package: 'Chem', name: 'getMpoProfileNames'})[0];
      return f ? ((await f.apply({})) as string[]) : [];
    })();
    if (names.length === 0) return; // no saved profiles on this stand — skip

    const columns = [
      DG.Column.fromFloat32Array('MW', new Float32Array([1, 2])),
      DG.Column.fromFloat32Array('LogP', new Float32Array([1, 2])),
    ];
    const harness = fakeCtx('', columns);
    const ed = factory(mappingParam(), harness.ctx);
    ed.setValue('');
    await until(() => ed.element.getAttribute('data-blocked') === 'true');

    harness.setProfile(names[0]);
    const built = await until(() =>
      ed.element.getAttribute('data-blocked') === null ||
      (ed.element.textContent ?? '').includes('scores no properties'), 20_000);
    expect(built, true, 'the editor reacted to the profile change on its own');
  });
});

category('Flow: func wrapper parity', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('a wrapped node renders the same parameters on the node and in the panel', async () => {
    for (const nqName of Object.keys(FUNC_WRAPPERS)) {
      const typeName = typeNameOf(nqName);
      if (!typeName) continue;
      const e = makeEditor();
      const panel = new PropertyPanel(e.flow);
      document.body.appendChild(panel.root);
      try {
        const node = await addNode(e.flow, typeName);
        const wrapper = funcWrapperOf(node.dgFunc!)!;
        const exposed = wrapper.inputs.map((i) => i.name);
        for (const name of exposed) {
          expect(name in node.inputs, true, `${nqName}: node has a socket for ${name}`);
          await until(() => !!e.container.querySelector(`[data-testid="${tid('socket-input', name)}"]`));
        }
        for (const raw of node.dgFunc!.inputs.map((p) => p.name)) {
          if (!exposed.includes(raw))
            expect(raw in node.inputs, false, `${nqName}: ${raw} is folded away, not a socket`);
        }
        panel.showNode(node);
        // a dataframe slot is connection-only — no data-param row
        for (const name of exposed) {
          if (name in node.inputValues)
            expect(!!panel.root.querySelector(`[data-param="${name}"]`), true, `${nqName}: panel row for ${name}`);
        }
        for (const raw of node.dgFunc!.inputs.map((p) => p.name)) {
          if (!exposed.includes(raw))
            expect(!!panel.root.querySelector(`[data-param="${raw}"]`), false,
              `${nqName}: no panel row for the folded-away ${raw}`);
        }
      } finally {
        panel.root.remove();
        destroyEditor(e);
      }
    }
  });

  test('effectiveFuncInputs returns the wrapper inputs, else the function\'s own', async () => {
    for (const nqName of Object.keys(FUNC_WRAPPERS)) {
      const typeName = typeNameOf(nqName);
      if (!typeName) continue;
      const func = getRegisteredFuncs().find((f) => f.nodeTypeName === typeName)!.func;
      expect(inputNames(func).join(','), FUNC_WRAPPERS[nqName].inputs.map((i) => i.name).join(','));
    }
    const plain = getRegisteredFuncs().find((f) => !funcWrapperOf(f.func) && f.func.inputs.length > 0);
    if (plain)
      expect(inputNames(plain.func).join(','), plain.func.inputs.map((p) => p.name).join(','), 'unwrapped: unchanged');
  });

  test('no wrapper invents a table input the function does not declare', async () => {
    for (const [nqName, wrapper] of Object.entries(FUNC_WRAPPERS)) {
      const typeName = typeNameOf(nqName);
      if (!typeName) continue;
      const func = getRegisteredFuncs().find((f) => f.nodeTypeName === typeName)!.func;
      const declared = new Set(func.inputs.map((p) => p.name));
      const invented = wrapper.inputs.filter((i) => !declared.has(i.name) && i.type === 'dataframe');
      const mapped = wrapper.mapInputs(Object.fromEntries(wrapper.inputs.map((i) => [i.name, `__${i.name}`])));
      for (const i of invented) {
        expect(Object.values(mapped).some((v) => String(v).includes(`__${i.name}`)), true,
          `${nqName}: ${i.name} is exposed but never reaches the call — use a real twin function instead`);
      }
    }
  });
});

category('Flow: hidden outputs', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('a hidden output keeps its slot but renders no socket row', async () => {
    const nqName = Object.keys(HIDDEN_FUNC_OUTPUTS)[0];
    const typeName = typeNameOf(nqName);
    if (!typeName) return;
    const e = makeEditor();
    try {
      const node = await addNode(e.flow, typeName);
      const hidden = Array.from(hiddenOutputsOf(node.dgFunc!));
      expect(hidden.length > 0, true, `${nqName} declares hidden outputs`);
      for (const key of hidden) {
        expect(key in node.outputs, true, `${key} slot still exists`);
        expect(node.hiddenOutputs.has(key), true, 'registry read onto the node');
      }
      // Wait for a row that IS expected, then assert the hidden ones are absent.
      await until(() => !!e.container.querySelector('[data-testid^="ff-socket-input"]'));
      for (const key of hidden) {
        expect(!!e.container.querySelector(`[data-testid="${tid('socket-output', key)}"]`), false,
          `no node row for hidden output ${key}`);
      }
    } finally {
      destroyEditor(e);
    }
  });

  test('hidden inputs and outputs never name the same key twice', async () => {
    for (const nqName of Object.keys(HIDDEN_FUNC_OUTPUTS)) {
      const ins = Object.keys(HIDDEN_FUNC_INPUTS[nqName] ?? {});
      for (const out of Object.keys(HIDDEN_FUNC_OUTPUTS[nqName]))
        expect(ins.includes(out), false, `${nqName}: ${out} is listed as both an input and an output`);
    }
  });
});

category('Flow: sketcher input', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('Sketcher Input is a string input carrying semType Molecule', async () => {
    const node = createNode('Inputs/Sketcher Input')!;
    expect(node !== null, true, 'registered');
    expect(node.dgOutputType, 'string', 'a sketched molecule IS a string — nothing downstream changes');
    expect(node.properties['semType'], 'Molecule');
  });

  test('inputValueProperty passes semType through to the built property', async () => {
    const sketcher = createNode('Inputs/Sketcher Input')!;
    expect(inputValueProperty(sketcher)?.semType, 'Molecule');
    const plain = createNode('Inputs/String Input')!;
    expect(inputValueProperty(plain)?.semType ?? '', '', 'a plain string input declares no semType');
    plain.properties['semType'] = 'Macromolecule';
    expect(inputValueProperty(plain)?.semType, 'Macromolecule', 'the qualifier is what drives it, not the node type');
  });

  test('a plain string drag still offers String Input first, not the sketcher', async () => {
    const producers = findNodeTypesProducingOutput('string');
    const inputs = producers.filter((p) => p.typeName.startsWith('Inputs/'));
    expect(inputs[0].typeName, 'Inputs/String Input',
      'a semType-specialized input never outranks the general one for a bare drag');
    expect(inputs.some((p) => p.typeName === 'Inputs/Sketcher Input'), true, 'but it is still offered');
  });

  test('the node body starts compact and expands into an in-node sketcher — no dialog', async () => {
    const node = createNode('Inputs/Sketcher Input')!;
    const ed = buildInputValueEditor(node, () => {}, {host: 'node'})!;
    document.body.appendChild(ed.root);
    try {
      expect(ed.root.classList.contains('ff-inline-sketcher'), true);
      const compact = ed.root.querySelector<HTMLElement>('[data-testid="ff-sketcher-compact"]')!;
      const full = ed.root.querySelector<HTMLElement>('[data-testid="ff-sketcher-full"]')!;
      expect(compact != null && full != null, true, 'both states exist');
      expect(full.style.display, 'none', 'starts compact');
      expect(compact.textContent, 'Sketch', 'an empty value shows the sketch invitation');
      expect(ed.sketcher == null, true, 'the sketcher is built lazily, on first expand');
      const dialogsBefore = DG.Dialog.getOpenDialogs().length;
      compact.click();
      expect(DG.Dialog.getOpenDialogs().length, dialogsBefore, 'expanding opens NO dialog');
      expect(ed.sketcher != null, true, 'the click built the inplace DG.chem.Sketcher');
      expect(full.style.display !== 'none', true, 'the sketcher state is shown');
      expect(compact.style.display, 'none', 'the compact preview folds away');
      expect(ed.root.classList.contains('ff-sketcher-expanded'), true);
      expect(ed.root.querySelector('[data-testid="ff-sketcher-done"]') != null, true, 'a Done control');
      (ed.root.querySelector<HTMLElement>('[data-testid="ff-sketcher-done"]'))!.click();
      expect(full.style.display, 'none', 'Done folds back to the compact preview');
      expect(compact.style.display !== 'none', true);
      // The panel Value row keeps the standard editor.
      const panelEd = buildInputValueEditor(node, () => {})!;
      expect(panelEd.root.classList.contains('ff-inline-sketcher'), false,
        'the panel Value row stays the standard molecule input');
      // The routing reads the semType qualifier, not the node type.
      const tagged = createNode('Inputs/String Input')!;
      tagged.properties['semType'] = 'Molecule';
      expect(hostsInlineSketcher(tagged), true, 'semType drives it');
      expect(hostsInlineSketcher(createNode('Inputs/Helm Input')!), false, 'Macromolecule routes to Helm, not chem');
    } finally {
      ed.root.remove();
    }
  });

  test('sketching writes SMILES into the configured value and reports the edit', async () => {
    const node = createNode('Inputs/Sketcher Input')!;
    let edits = 0;
    const ed = buildInputValueEditor(node, () => edits++, {host: 'node'})!;
    document.body.appendChild(ed.root);
    try {
      ed.root.querySelector<HTMLElement>('[data-testid="ff-sketcher-compact"]')!.click();
      // A sketch gesture = a real interaction followed by the sketcher's change event.
      // setSmiles, not setMolecule — the latter routes through the synchronous
      // Chem:isSmarts, which throws when the harness hasn't loaded Chem yet.
      ed.root.dispatchEvent(new PointerEvent('pointerdown', {bubbles: true}));
      ed.sketcher!.setSmiles('c1ccccc1');
      ed.sketcher!.onChanged.next(null);
      expect(String(node.properties['defaultValue'] ?? '') !== '', true, 'the SMILES landed in the value');
      expect(edits, 1, 'exactly one reported edit');
      const r = resolveInputValue(node);
      expect(r.ok, true, 'the run can feed the sketched molecule to the prepared call');
      expect(String(r.value).includes('\n'), false, 'SMILES, never a multiline molfile (header emission)');
      // The same molecule again is not another edit.
      ed.sketcher!.onChanged.next(null);
      expect(edits, 1, 'an unchanged value never re-reports');
    } finally {
      ed.root.remove();
    }
  });

  test('programmatic loads and syncs never report an edit', async () => {
    const node = createNode('Inputs/Sketcher Input')!;
    node.properties['defaultValue'] = 'CCO';
    let edits = 0;
    const ed = buildInputValueEditor(node, () => edits++, {host: 'node'})!;
    document.body.appendChild(ed.root);
    try {
      ed.root.querySelector<HTMLElement>('[data-testid="ff-sketcher-compact"]')!.click();
      // The sketcher echoes programmatic sets through onChanged (async, canonicalized) —
      // untouched, those echoes must stay silent or merely loading a flow marks it dirty.
      ed.sketcher!.onChanged.next(null);
      expect(edits, 0, 'the load echo is not an edit');
      node.properties['defaultValue'] = 'CCC';
      ed.sync(); // the panel wrote a new value — programmatic too
      ed.sketcher!.onChanged.next(null);
      expect(edits, 0, 'the sync echo is not an edit');
      expect(ed.root.classList.contains('ff-value-missing'), false, 'a configured value is not flagged');
      node.properties['defaultValue'] = '';
      ed.sync();
      expect(ed.root.classList.contains('ff-value-missing'), true, 'an empty required value gets the amber cue');
    } finally {
      ed.root.remove();
    }
  });

  test('expanding snaps to native zoom, shows the whole 500×500 sketcher; a user zoom folds it', async () => {
    // The sketcher only exists at zoom 1 — its chrome is not built to be scaled.
    const e = makeEditor();
    try {
      // Bottom-right-ish so the visibility pan actually has work to do.
      const node = await addNode(e.flow, 'Inputs/Sketcher Input', 700, 400);
      const compactSel = `.ff-node[data-node-id="${node.id}"] [data-testid="ff-sketcher-compact"]`;
      expect(await until(() => e.container.querySelector(compactSel) != null, 10000), true,
        'the compact preview renders in the card');
      expect(estimateNodeWidth(node) <= 280, true, 'compact never lifts the card cap');
      e.flow.setZoom(0.5); // the user had zoomed out before clicking
      expect(await until(() => Math.abs(e.flow.getZoom() - 0.5) < 0.001, 2000), true, 'zoomed out');
      e.container.querySelector<HTMLElement>(compactSel)!.click();
      expect(await until(() => Math.abs(e.flow.getZoom() - 1) < 0.001, 3000), true,
        'expanding snaps the canvas to native zoom');
      const full = e.container.querySelector<HTMLElement>(
        `.ff-node[data-node-id="${node.id}"] [data-testid="ff-sketcher-full"]`)!;
      expect(full != null && full.style.display !== 'none', true, 'expanded in place');
      const r = full.getBoundingClientRect();
      expect(Math.abs(r.width - 500) <= 3 && Math.abs(r.height - 500) <= 3, true,
        `a fixed 500×500 box (was ${r.width}×${r.height})`);
      // The visibility pan is async — poll until the box settles inside the viewport.
      const inViewport = (): boolean => {
        const b = full.getBoundingClientRect();
        const c = e.container.querySelector<HTMLElement>('.ff-canvas')!.getBoundingClientRect();
        return b.left >= c.left - 1 && b.top >= c.top - 1 &&
          b.right <= c.right + 1 && b.bottom <= c.bottom + 1;
      };
      const ok = await until(inViewport, 3000);
      const b = full.getBoundingClientRect();
      const c = e.container.querySelector<HTMLElement>('.ff-canvas')!.getBoundingClientRect();
      expect(ok, true, 'the sketcher is panned fully into the viewport — box ' +
        `${JSON.stringify({l: b.left, t: b.top, r: b.right, b: b.bottom, display: full.style.display})} vs canvas ` +
        `${JSON.stringify({l: c.left, t: c.top, r: c.right, b: c.bottom})} zoom ${e.flow.getZoom()}`);
      // A user zoom hides the sketcher — back to the compact preview.
      e.flow.setZoom(1.5);
      expect(await until(() => full.style.display === 'none', 3000), true,
        'zooming away folds the sketcher back');
      const compact = e.container.querySelector<HTMLElement>(compactSel)!;
      expect(compact.style.display !== 'none', true, 'the compact preview is back');
    } finally {
      destroyEditor(e);
    }
  }, {timeout: 30000});
});

category('Flow: chem nodes', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  const TABLE_COLUMN_TWINS = [
    {nqName: 'Chem:filterBySubstructure', column: 'molecules', extra: 'substructure'},
    {nqName: 'Chem:similarityTo', column: 'molecules', extra: 'query'},
    {nqName: 'Chem:diverseSubset', column: 'molecules', extra: 'limit'},
    {nqName: 'Chem:applyReaction', column: 'molecules', extra: 'reaction'},
    {nqName: 'Chem:toSdf', column: 'molecules', extra: null},
  ];

  test('each chem twin takes a table plus a Molecule column, both required', async () => {
    for (const spec of TABLE_COLUMN_TWINS) {
      const typeName = typeNameOf(spec.nqName);
      if (!typeName) continue;
      const e = makeEditor();
      try {
        const node = await addNode(e.flow, typeName);
        const func = node.dgFunc!;
        const inputs = effectiveFuncInputs(func);
        expect(String(inputs[0].propertyType), 'dataframe', `${spec.nqName}: leads with a table`);
        const col = inputs.find((p) => p.name === spec.column)!;
        expect(col !== undefined, true, `${spec.nqName}: declares ${spec.column}`);
        expect(String(col.propertyType), 'column');
        expect(col.semType, 'Molecule', `${spec.nqName}: ${spec.column} is semType-filtered`);
        const tables = node.properties['columnTables'] as Record<string, string> | undefined;
        expect(tables?.[spec.column], inputs[0].name, `${spec.nqName}: ${spec.column} resolves against the table`);
        expect(node.requiredInputs.includes(inputs[0].name), true, `${spec.nqName}: table required`);
        expect(node.requiredInputs.includes(spec.column), true, `${spec.nqName}: ${spec.column} required`);
        if (spec.extra)
          expect(spec.extra in node.inputs, true, `${spec.nqName}: ${spec.extra} slot exists`);
      } finally {
        destroyEditor(e);
      }
    }
  });

  test('the chem twin panel offers a column picker for its molecule column', async () => {
    const typeName = typeNameOf('Chem:filterBySubstructure');
    if (!typeName) return;
    const e = makeEditor();
    const panel = new PropertyPanel(e.flow);
    panel.onPickColumns = () => {}; // the icon is gated on a handler being wired
    document.body.appendChild(panel.root);
    try {
      const node = await addNode(e.flow, typeName);
      panel.showNode(node);
      expect(!!panel.root.querySelector('[data-param="molecules"]'), true, 'the column row renders');
      expect(!!panel.root.querySelector(`[data-testid="${tid('prop-pick-columns', 'molecules')}"]`), true,
        'and carries the picker icon — a bare "connected only" label is the failure mode');
      expect(!!panel.root.querySelector('[data-param="substructure"]'), true, 'the query row renders');
    } finally {
      panel.root.remove();
      destroyEditor(e);
    }
  });

  test('the superseded column-only chem functions are off the allowlist', async () => {
    for (const gone of ['Chem:getSimilarities', 'Chem:getDiversities', 'Chem:searchSubstructure',
      'Chem:getInchis', 'Chem:getInchiKeys', 'Chem:getMorganFingerprints'])
      expect(INCLUDED_FUNC_NQNAMES.has(gone), false, `${gone} superseded`);
    for (const present of ['Chem:filterBySubstructure', 'Chem:similarityTo', 'Chem:diverseSubset',
      'Chem:addInchisTopMenu', 'Chem:addInchisKeysTopMenu'])
      expect(INCLUDED_FUNC_NQNAMES.has(present), true, `${present} is the replacement`);
  });

  test('Filter by Substructure really filters (needs Chem deployed)', async () => {
    const func = DG.Func.find({package: 'Chem', name: 'filterBySubstructure'})[0];
    if (!func) return;
    const molecules = DG.Column.fromStrings('smiles', ['c1ccccc1', 'CCO', 'c1ccccc1C', 'CCN']);
    molecules.semType = DG.SEMTYPE.MOLECULE;
    const table = DG.DataFrame.fromColumns([molecules, DG.Column.fromInt32Array('id', new Int32Array([1, 2, 3, 4]))]);
    const res: DG.DataFrame = await func.apply({table, molecules, substructure: 'c1ccccc1'});
    expect(res instanceof DG.DataFrame, true, 'returns a dataframe');
    expect(res.rowCount, 2, 'the two aromatics match, the two aliphatics do not');
    expect(res.columns.names().includes('id'), true, 'and the other columns come along');
    expect(table.rowCount, 4, 'the input table is not mutated');
  });

  test('Similarity To adds a score column to the table (needs Chem deployed)', async () => {
    const func = DG.Func.find({package: 'Chem', name: 'similarityTo'})[0];
    if (!func) return;
    const molecules = DG.Column.fromStrings('smiles', ['c1ccccc1', 'CCO', 'c1ccccc1C']);
    molecules.semType = DG.SEMTYPE.MOLECULE;
    const table = DG.DataFrame.fromColumns([molecules]);
    const col: DG.Column = await func.apply({table, molecules, query: 'c1ccccc1'});
    expect(col instanceof DG.Column, true, 'returns a column');
    expect(col.length, 3, 'one score per row');
    // exact column count not asserted — the fingerprint cache attaches a column of its own
    expect(table.col(col.name) !== null, true, 'the score column is in the table under the name it reports');
    expect(col.stats.max <= 1 && col.stats.min >= 0, true, 'Tanimoto scores are in 0..1');
  });
});
