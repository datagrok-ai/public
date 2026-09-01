import * as DG from 'datagrok-api/dg';
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {
  registerBuiltinNodes, registerAllFunctions, createNode, ensureFuncNodeType, getRegisteredTypeNames,
  getRegisteredFuncs, findNodeTypesAcceptingInput, findNodeTypesProducingOutput, funcCategory,
  candidateMatchesQuery, prioritizeCandidates, CompatibleNodeType,
} from '../rete/node-factory';
import {FUNC_CATEGORIES} from '../panel/function-browser';
import {FuncNode} from '../rete/nodes/func-node';
import {getParamDescription, getParamDisplayName} from '../utils/dart-proxy-utils';
import {pastelize, FUNC_NAME_COLORS} from '../types/type-map';
import {makeEditor, destroyEditor, addNode, until} from './test-utils';

category('Flow: node-factory', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('createNode builds known built-in types', async () => {
    const table = createNode('Inputs/Table Input');
    expect(table !== null, true);
    expect(table!.dgNodeType, 'input');
    expect(table!.dgTypeName, 'Inputs/Table Input');
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

  test('a leading (table, column) input pair is required even when annotated nullable', async () => {
    const func = DG.Func.find({package: 'Chem', name: 'bitbirchClusteringTopMenu'})[0];
    if (!func) return; // Chem not on this stand
    const typeName = ensureFuncNodeType(func);
    const node = createNode(typeName)!;
    expect(node.requiredInputs.includes('table'), true, 'nullable leading table is forced required');
    expect(node.requiredInputs.includes('molecules'), true, 'nullable leading column is forced required');
    expect(node.requiredInputs.includes('threshold'), false, 'defaulted param stays optional');

    const open = getRegisteredFuncs().find((f) => f.func.name === 'OpenFile');
    if (open) {
      const of = createNode(open.nodeTypeName)!;
      expect(of.requiredInputs.includes('fullPath'), true, 'genuinely required stays required');
      expect(of.requiredInputs.includes('sheetName'), false, 'heuristic off without the pair');
    }
  });

  test('built-in registry includes the expected core types', async () => {
    const names = new Set(getRegisteredTypeNames());
    for (const t of ['Inputs/Table Input', 'Outputs/Table Output', 'Outputs/Value Output',
      'Constants/String', 'Constants/Int', 'Constants/Boolean'])
      expect(names.has(t), true, `registry has ${t}`);
  });

  test('SetVar and GetVar are always registered — saved flows depend on them', async () => {
    if (DG.Func.find({name: 'SetVar'}).length === 0) {
      expect(true, true); // no live backend — nothing to check
      return;
    }
    for (const name of ['setvar', 'getvar']) {
      const info = getRegisteredFuncs().find((f) => (f.func.name ?? '').toLowerCase() === name);
      expect(!!info, true, `${name} is in the registry after registerAllFunctions()`);
      const node = createNode(info!.nodeTypeName);
      expect(!!node, true, `${name} node type instantiates`);
      expect(node!.dgNodeType, 'func');
    }
  });

  test('ensureFuncNodeType is idempotent for the same function', async () => {
    const func = DG.Func.find({name: 'AddNewColumn'})[0] ?? DG.Func.find({})[0];
    expect(func != null, true);
    const a = ensureFuncNodeType(func);
    const b = ensureFuncNodeType(func);
    expect(a, b);
    const node = createNode(a);
    expect(node instanceof FuncNode, true);
    expect(node!.dgFunc?.name, func.name);
  });

  test('FuncNode builds a pass-through output per input', async () => {
    const func = DG.Func.find({name: 'AddNewColumn'})[0];
    if (!func) return;
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
    const setVarRed = FUNC_NAME_COLORS['setvar'].color;
    const setVar = DG.Func.find({name: 'SetVar'})[0];
    if (!setVar) return;
    const node = new FuncNode(setVar);
    expect((node as unknown as {color?: string}).color, setVarRed, 'SetVar title is red');

    const other = DG.Func.find({name: 'AddNewColumn'})[0];
    if (other)
      expect((new FuncNode(other) as unknown as {color?: string}).color !== setVarRed, true);
  });

  test('title bar renders the pastel of the identity color, not the vivid hue', async () => {
    const e = makeEditor();
    try {
      const a = await addNode(e.flow, 'Utilities/Info', 0, 0);
      const query = (): HTMLElement | null => e.container.querySelector<HTMLElement>(
        `.ff-node[data-node-id="${a.id}"] .ff-node-title`);
      await until(() => query() != null);
      const title = query()!;
      // Normalize both sides through the browser's own color parsing.
      const norm = (c: string): string => {
        const d = document.createElement('div');
        d.style.background = c;
        return d.style.background;
      };
      const identity = (a as unknown as {color: string}).color;
      expect(norm(identity) !== '', true, 'node keeps a vivid identity color');
      expect(title.style.background, norm(pastelize(identity)), 'title bar is the pastel');
      expect(title.style.background !== norm(identity), true, 'vivid hue is not painted directly');
    } finally {
      destroyEditor(e);
    }
  });

  test('sockets render Column Manager-style type letters; order squares stay bare', async () => {
    const e = makeEditor();
    try {
      const a = await addNode(e.flow, 'Inputs/Table Input', 0, 0);
      const sel = `.ff-node[data-node-id="${a.id}"]`;
      await until(() => e.container.querySelector(`${sel} [data-testid="ff-socket-dataframe"]`) != null);
      const chip = e.container.querySelector<HTMLElement>(`${sel} [data-testid="ff-socket-dataframe"]`)!;
      expect(chip.textContent, 't', 'a dataframe socket shows the "t" (table) letter');
      expect(chip.style.getPropertyValue('--socket-color') !== '', true, 'chip carries its type color');
      const order = e.container.querySelector<HTMLElement>(`${sel} .ff-exec-out .ff-socket`)!;
      expect(order.textContent, '', 'order squares carry no letter');
    } finally {
      destroyEditor(e);
    }
  });

  test('suggestion menu shows friendly names with "what it does" categories', async () => {
    const candidates = findNodeTypesAcceptingInput('dataframe');
    expect(candidates.length > 0, true, 'a table output has compatible next steps');

    const uncategorized = candidates.filter((c) => c.label.includes('Uncategorized'));
    expect(uncategorized.length, 0,
      `"Uncategorized" leaked into labels: ${uncategorized.slice(0, 3).map((c) => c.label).join(' | ')}`);

    const anc = getRegisteredFuncs().find((f) => f.func.name === 'AddNewColumn');
    if (anc) {
      const item = candidates.find((c) => c.typeName === anc.nodeTypeName);
      expect(!!item, true, 'AddNewColumn is offered for a table output');
      expect(item!.label, `${anc.name}  (${funcCategory(anc)})`, 'friendly name + what-it-does category');
      expect(funcCategory(anc), 'Column Operations', 'AddNewColumn is a column operation');
    }

    // Only the TRAILING paren group is the category — a friendly name may itself
    // carry parentheses ("Chemical Properties (OCL)"), so a greedy match is wrong.
    const catRe = /\(([^()]+)\)$/;
    for (const c of candidates.filter((x) => !x.isBuiltin).slice(0, 50)) {
      const m = catRe.exec(c.label);
      expect(!!m, true, `label has a category: ${c.label}`);
      expect((FUNC_CATEGORIES as readonly string[]).includes(m![1]), true, `known category in ${c.label}`);
    }

    const builtin = candidates.find((c) => c.typeName === 'Outputs/Table Output');
    if (builtin) expect(builtin.label, 'Table Output');

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

    const plain = findNodeTypesAcceptingInput('dataframe');
    expect(plain[0].typeName, 'Outputs/Value Output', 'Value Output leads');

    const tableOut = idxOf(plain, (c) => c.typeName === 'Outputs/Table Output');
    const log = idxOf(plain, (c) => c.typeName === 'Utilities/Log');
    if (tableOut !== -1 && log !== -1)
      expect(tableOut < log, true, 'exact dataframe consumer before a dynamic catch-all');

    if (getRegisteredFuncs().some((f) => f.packageName === 'Chem')) {
      const fromChem = findNodeTypesAcceptingInput('dataframe', {sourcePackageName: 'Chem'});
      const lastChem = (fromChem as Array<{label: string}>).map(isChem).lastIndexOf(true);
      const ancFrom = idxOf(fromChem, byFunc('AddNewColumn'));
      expect(lastChem !== -1, true, 'chem candidates exist for a table drag');
      if (ancFrom !== -1)
        expect(lastChem < ancFrom, true, 'all chem funcs precede the common core funcs');

      const viaGraph = findNodeTypesAcceptingInput('dataframe', {graphPackageNames: ['Chem']});
      expect(isChem(viaGraph[1] as {label: string}), true, 'canvas domain boosts chem right after Value Output');

      const ancPlain = idxOf(plain, byFunc('AddNewColumn'));
      const firstChemPlain = idxOf(plain, isChem as (c: {typeName: string; label: string}) => boolean);
      if (ancPlain !== -1 && firstChemPlain !== -1)
        expect(ancPlain < firstChemPlain, true, 'without context, chem is not boosted');
    }

    const ancIdx = idxOf(plain, byFunc('AddNewColumn'));
    const aggIdx = idxOf(plain, byFunc('Aggregate'));
    if (ancIdx !== -1 && aggIdx !== -1) {
      expect(ancIdx < aggIdx, true, 'baseline: Add New Column sorts before Aggregate');
      const used = findNodeTypesAcceptingInput('dataframe', {graphFuncNames: ['Aggregate']});
      expect(idxOf(used, byFunc('Aggregate')) < idxOf(used, byFunc('AddNewColumn')), true,
        'a func already on the canvas floats above its tier peers');
    }
  });

  test('suggestion menu search covers descriptions, tags, and package like the toolbox', async () => {
    const c: CompatibleNodeType = {
      typeName: 'DG Functions/Uncategorized/Flow:deleteColumns',
      label: 'Delete Columns  (Transform Tables)', isBuiltin: false,
      searchText: 'delete columns  (transform tables) flow removes the given columns, as a new table',
    };
    expect(candidateMatchesQuery(c, ''), true, 'empty query matches everything');
    expect(candidateMatchesQuery(c, 'remove'), true, 'description word matches');
    expect(candidateMatchesQuery(c, 'removes the given'), true, 'description phrase matches');
    expect(candidateMatchesQuery(c, 'deletecolumns'), true, 'whitespace-insensitive name match');
    expect(candidateMatchesQuery(c, 'Flow'), true, 'package matches, case-insensitive');
    expect(candidateMatchesQuery(c, 'zzz-nothing'), false, 'non-matching query rejected');
    expect(candidateMatchesQuery({typeName: 'Utilities/Log', label: 'Log', isBuiltin: true}, 'log'),
      true, 'haystack-less candidate falls back to label/typeName');

    const candidates = findNodeTypesAcceptingInput('dataframe');
    const tableOut = candidates.find((x) => x.typeName === 'Outputs/Table Output');
    expect(!!tableOut?.searchText, true, 'built-in candidates carry a search haystack');
    expect(candidateMatchesQuery(tableOut!, 'marks a dataframe'), true,
      'built-in description searchable (same text as the toolbox tooltip)');

    const withDesc = getRegisteredFuncs().find((f) => {
      try {
        return String(f.func.description || '').length > 10 &&
          candidates.some((x) => x.typeName === f.nodeTypeName);
      } catch {
        return false;
      }
    });
    if (withDesc) {
      const item = candidates.find((x) => x.typeName === withDesc.nodeTypeName)!;
      const fragment = String(withDesc.func.description).toLowerCase().slice(0, 12);
      expect(candidateMatchesQuery(item, fragment), true,
        `func description searchable in the menu: "${fragment}" → ${item.label}`);
    }

    const producers = findNodeTypesProducingOutput('dataframe');
    const tableIn = producers.find((x) => x.typeName === 'Inputs/Table Input');
    expect(candidateMatchesQuery(tableIn!, 'dataframe input parameter'), true,
      'reverse-menu built-in description searchable');
  });

  test('prioritizeCandidates floats the suggestion-engine picks, in engine order, with reasons', async () => {
    const cand = (typeName: string): CompatibleNodeType => ({typeName, label: typeName, isBuiltin: true});
    const candidates = [cand('A'), cand('B'), cand('C'), cand('D')];

    const merged = prioritizeCandidates(candidates, [
      {typeName: 'C', reason: 'Molecule column "smiles"'},
      {typeName: 'B', reason: 'Table from "Open File"', prefill: {molecules: 'smiles'}},
      {typeName: 'X', reason: 'not a menu candidate — dropped'},
    ]);
    expect(merged.map((x) => x.typeName).join(','), 'C,B,A,D',
      'engine picks lead in engine order, the rest keep their order');
    expect(merged[0].reason, 'Molecule column "smiles"', 'reason attached');
    expect(merged[2].reason == null, true, 'non-suggested items carry no reason');
    expect(candidates[2].reason == null, true, 'input list not mutated');

    expect(prioritizeCandidates(candidates, []), candidates, 'no suggestions → same list');
  });

  test('reverse suggestions: producers of a type, real outputs before passthrough threaders', async () => {
    const idxOf = (list: {typeName: string}[], pred: (c: {typeName: string; label: string}) => boolean): number =>
      (list as Array<{typeName: string; label: string}>).findIndex(pred);
    const byFunc = (name: string) => (c: {typeName: string}): boolean =>
      (c.typeName.split('/').pop() ?? '').split(':').pop() === name;

    const tables = findNodeTypesProducingOutput('dataframe');
    expect(tables.length > 0, true, 'table producers exist');
    expect(tables[0].typeName, 'Inputs/Table Input', 'the matching Input node leads');

    const openFile = idxOf(tables, byFunc('OpenFile'));
    if (openFile !== -1) {
      const threader = tables.findIndex((c) => c.realOutput === false);
      if (threader !== -1)
        expect(openFile < threader, true, 'a real producer precedes passthrough-only threaders');
    }

    const agg = idxOf(tables, byFunc('Aggregate'));
    const anc = idxOf(tables, byFunc('AddNewColumn'));
    if (agg !== -1 && anc !== -1) {
      expect(agg < anc, true, 'real dataframe output beats a passthrough-only match');
      const ancItem = tables[anc] as {realOutput?: boolean};
      expect(ancItem.realOutput, false, 'AddNewColumn matches via passthrough only');
    }

    const strings = findNodeTypesProducingOutput('string');
    expect(strings[0].typeName, 'Inputs/String Input', 'String Input leads a string drag');

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
