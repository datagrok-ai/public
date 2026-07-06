import * as DG from 'datagrok-api/dg';
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {
  registerBuiltinNodes, registerAllFunctions, createNode, ensureFuncNodeType, getRegisteredTypeNames,
  getRegisteredFuncs, findNodeTypesAcceptingInput, findNodeTypesProducingOutput, funcCategory,
} from '../rete/node-factory';
import {FUNC_CATEGORIES} from '../panel/function-browser';
import {FuncNode} from '../rete/nodes/func-node';
import {getParamDescription, getParamDisplayName} from '../utils/dart-proxy-utils';

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

  test('FuncNode records the source package (shown in the context panel)', async () => {
    const info = getRegisteredFuncs().find((f) => f.packageName === 'Chem') ??
      getRegisteredFuncs().find((f) => !!f.packageName);
    if (!info) return; // no package funcs on this stand — skip
    const node = new FuncNode(info.func);
    expect(node.dgPackageName, info.packageName, 'package name captured on the node');
  });

  test('input slot label shows the caption, key stays the property name', async () => {
    // Find a registered func whose input declares a caption distinct from its
    // name (e.g. Aggregate/PCA/dbScan) — the slot label should be that caption
    // while the slot key stays the property name.
    let found: {func: DG.Func; name: string; caption: string} | null = null;
    for (const info of getRegisteredFuncs()) {
      for (const p of info.func.inputs) {
        const cap = getParamDisplayName(p);
        if (cap && cap !== p.name) {found = {func: info.func, name: p.name, caption: cap}; break;}
      }
      if (found) break;
    }
    if (!found) return; // no caption-bearing inputs on this stand — skip
    const node = new FuncNode(found.func);
    expect(found.name in node.inputs, true, 'slot key is the property name (identity unchanged)');
    expect((node.inputs[found.name] as {label?: string}).label, found.caption, 'slot label is the caption');
  });

  test('FuncNode captures param descriptions for socket/panel tooltips', async () => {
    // Find any registered func whose inputs OR outputs declare a description
    // (via @grok.decorators.param) — validates the Dart-proxy read end to end.
    let found: {func: DG.Func; param: string} | null = null;
    for (const info of getRegisteredFuncs()) {
      for (const p of [...info.func.inputs, ...info.func.outputs])
        if (getParamDescription(p)) {found = {func: info.func, param: p.name}; break;}
      if (found) break;
    }
    if (!found) return; // no described params on this stand — skip
    const node = new FuncNode(found.func);
    const all = {...node.inputDescriptions, ...node.outputDescriptions};
    expect(found.param in all, true, `description captured for ${found.func.name}.${found.param}`);
    expect(all[found.param].length > 0, true, 'description is non-empty');
  });

  test('per-function color override pins SetVar to red', async () => {
    const setVar = DG.Func.find({name: 'SetVar'})[0];
    if (!setVar) return;
    const node = new FuncNode(setVar);
    expect((node as unknown as {color?: string}).color, '#EF5350', 'SetVar title is red');

    // A function without an override falls back to role/default coloring.
    const other = DG.Func.find({name: 'AddNewColumn'})[0];
    if (other)
      expect((new FuncNode(other) as unknown as {color?: string}).color !== '#EF5350', true);
  });

  test('suggestion menu shows friendly names with "what it does" categories', async () => {
    // Dragging a table output to empty canvas lists compatible nodes. DG funcs
    // must read like the toolbox — friendly name + task category — not the raw
    // `funcName (Uncategorized)` baked into the typeName.
    const candidates = findNodeTypesAcceptingInput('dataframe');
    expect(candidates.length > 0, true, 'a table output has compatible next steps');

    // No DG function is labeled with the role-segment fallback.
    const uncategorized = candidates.filter((c) => c.label.includes('Uncategorized'));
    expect(uncategorized.length, 0,
      `"Uncategorized" leaked into labels: ${uncategorized.slice(0, 3).map((c) => c.label).join(' | ')}`);

    // A known catalog func shows its registry display name and task category.
    const anc = getRegisteredFuncs().find((f) => f.func.name === 'AddNewColumn');
    if (anc) {
      const item = candidates.find((c) => c.typeName === anc.nodeTypeName);
      expect(!!item, true, 'AddNewColumn is offered for a table output');
      expect(item!.label, `${anc.name}  (${funcCategory(anc)})`, 'friendly name + what-it-does category');
      expect(funcCategory(anc), 'Column Operations', 'AddNewColumn is a column operation');
    }

    // Every DG-func candidate carries a known category in parentheses.
    const catRe = /\((.+)\)$/;
    for (const c of candidates.filter((x) => !x.isBuiltin).slice(0, 50)) {
      const m = catRe.exec(c.label);
      expect(!!m, true, `label has a category: ${c.label}`);
      expect((FUNC_CATEGORIES as readonly string[]).includes(m![1]), true, `known category in ${c.label}`);
    }

    // Built-ins keep their plain label.
    const builtin = candidates.find((c) => c.typeName === 'Outputs/Table Output');
    if (builtin) expect(builtin.label, 'Table Output');

    // A Chem operation lands under its domain, same as the toolbox.
    const chemOp = candidates.find((c) => {
      if (c.isBuiltin) return false;
      const info = getRegisteredFuncs().find((f) => f.nodeTypeName === c.typeName);
      return info?.packageName === 'Chem';
    });
    if (chemOp) expect(chemOp.label.endsWith('(Cheminformatics)'), true, `chem op labeled by domain: ${chemOp.label}`);
  });

  test('suggestion ranking: the science in play leads, exact types beat wildcards, used funcs float', async () => {
    const idxOf = (list: {typeName: string}[], pred: (c: {typeName: string; label: string}) => boolean): number =>
      (list as Array<{typeName: string; label: string}>).findIndex(pred);
    const isChem = (c: {label: string}): boolean => c.label.endsWith('(Cheminformatics)');
    const byFunc = (name: string) => (c: {typeName: string}): boolean =>
      (c.typeName.split('/').pop() ?? '').split(':').pop() === name;

    // Value Output is always first, context or not.
    const plain = findNodeTypesAcceptingInput('dataframe');
    expect(plain[0].typeName, 'Outputs/Value Output', 'Value Output leads');

    // Exact type match beats a wildcard acceptor in the same tier: for a table
    // drag, Table Output (dataframe input) precedes Log (dynamic input) even
    // though "Log" sorts first alphabetically.
    const tableOut = idxOf(plain, (c) => c.typeName === 'Outputs/Table Output');
    const log = idxOf(plain, (c) => c.typeName === 'Utilities/Log');
    if (tableOut !== -1 && log !== -1)
      expect(tableOut < log, true, 'exact dataframe consumer before a dynamic catch-all');

    if (getRegisteredFuncs().some((f) => f.packageName === 'Chem')) {
      // Dragging out of a Chem node → every Cheminformatics function precedes
      // the common core funcs (and thus all other DG funcs).
      const fromChem = findNodeTypesAcceptingInput('dataframe', {sourcePackageName: 'Chem'});
      const lastChem = (fromChem as Array<{label: string}>).map(isChem).lastIndexOf(true);
      const ancFrom = idxOf(fromChem, byFunc('AddNewColumn'));
      expect(lastChem !== -1, true, 'chem candidates exist for a table drag');
      if (ancFrom !== -1)
        expect(lastChem < ancFrom, true, 'all chem funcs precede the common core funcs');

      // A domain-less source (OpenFile, a utility) falls back to the science
      // already on the canvas.
      const viaGraph = findNodeTypesAcceptingInput('dataframe', {graphPackageNames: ['Chem']});
      expect(isChem(viaGraph[1] as {label: string}), true, 'canvas domain boosts chem right after Value Output');

      // No context → no domain boost: the common core funcs stay ahead of chem.
      const ancPlain = idxOf(plain, byFunc('AddNewColumn'));
      const firstChemPlain = idxOf(plain, isChem as (c: {typeName: string; label: string}) => boolean);
      if (ancPlain !== -1 && firstChemPlain !== -1)
        expect(ancPlain < firstChemPlain, true, 'without context, chem is not boosted');
    }

    // A function already used on the canvas floats within its tier: Aggregate
    // jumps ahead of Add New Column (which otherwise wins alphabetically).
    const ancIdx = idxOf(plain, byFunc('AddNewColumn'));
    const aggIdx = idxOf(plain, byFunc('Aggregate'));
    if (ancIdx !== -1 && aggIdx !== -1) {
      expect(ancIdx < aggIdx, true, 'baseline: Add New Column sorts before Aggregate');
      const used = findNodeTypesAcceptingInput('dataframe', {graphFuncNames: ['Aggregate']});
      expect(idxOf(used, byFunc('Aggregate')) < idxOf(used, byFunc('AddNewColumn')), true,
        'a func already on the canvas floats above its tier peers');
    }
  });

  test('reverse suggestions: producers of a type, real outputs before passthrough threaders', async () => {
    const idxOf = (list: {typeName: string}[], pred: (c: {typeName: string; label: string}) => boolean): number =>
      (list as Array<{typeName: string; label: string}>).findIndex(pred);
    const byFunc = (name: string) => (c: {typeName: string}): boolean =>
      (c.typeName.split('/').pop() ?? '').split(':').pop() === name;

    // "What produces a table?" — the matching Input node leads (the universal
    // "make this a script parameter" producer).
    const tables = findNodeTypesProducingOutput('dataframe');
    expect(tables.length > 0, true, 'table producers exist');
    expect(tables[0].typeName, 'Inputs/Table Input', 'the matching Input node leads');

    // A Data Sources func (OpenFile — real dataframe output) is boosted.
    const openFile = idxOf(tables, byFunc('OpenFile'));
    if (openFile !== -1) {
      // It precedes any passthrough-only threader (e.g. a func that merely
      // threads a table through, like a column-outputting mutator).
      const threader = tables.findIndex((c) => c.realOutput === false);
      if (threader !== -1)
        expect(openFile < threader, true, 'a real producer precedes passthrough-only threaders');
    }

    // Real-over-passthrough within a tier: Aggregate (real dataframe output)
    // beats Add New Column (only its table passthrough matches a table drag) —
    // both are COMMON funcs, and alphabetics would put ANC first.
    const agg = idxOf(tables, byFunc('Aggregate'));
    const anc = idxOf(tables, byFunc('AddNewColumn'));
    if (agg !== -1 && anc !== -1) {
      expect(agg < anc, true, 'real dataframe output beats a passthrough-only match');
      const ancItem = tables[anc] as {realOutput?: boolean};
      expect(ancItem.realOutput, false, 'AddNewColumn matches via passthrough only');
    }

    // A string drag leads with String Input; a passthrough-only threader for
    // strings (any func with a string input) is offered too, ranked below.
    const strings = findNodeTypesProducingOutput('string');
    expect(strings[0].typeName, 'Inputs/String Input', 'String Input leads a string drag');

    // The domain boost applies in the reverse direction too.
    if (getRegisteredFuncs().some((f) => f.packageName === 'Chem')) {
      const fromChem = findNodeTypesProducingOutput('dataframe', {sourcePackageName: 'Chem'});
      const isChem = (c: {label: string}): boolean => c.label.endsWith('(Cheminformatics)');
      const firstChem = idxOf(fromChem, isChem as (c: {typeName: string; label: string}) => boolean);
      const ancChem = idxOf(fromChem, byFunc('Aggregate'));
      if (firstChem !== -1 && ancChem !== -1)
        expect(firstChem < ancChem, true, 'chem producers lead when dragging from a chem node');
    }
  });
});
