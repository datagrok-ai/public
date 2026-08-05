import * as DG from 'datagrok-api/dg';
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes, registerAllFunctions, ensureFuncNodeType} from '../rete/node-factory';
import {isExecKey} from '../rete/scheme';
import {buildCreationScriptGraph, applyGraphToEditor} from '../import/creation-script-importer';
import {
  emitCreationScript, emitCreationScriptsForTables,
  CreationScriptResult, PerTableResult,
} from '../compiler/creation-script-emitter';
import {makeEditor, destroyEditor, addNode} from './test-utils';

function chemAvailable(): boolean {
  try {
    return DG.Func.find({name: 'addChemPropertiesColumns'}).length > 0;
  } catch {
    return false;
  }
}

/** Import a creation script into a live editor and re-emit it. */
async function roundTrip(script: string): Promise<CreationScriptResult> {
  const e = makeEditor();
  try {
    const g = buildCreationScriptGraph(script);
    await applyGraphToEditor(g, e.flow);
    return emitCreationScript(e.flow);
  } finally {
    destroyEditor(e);
  }
}

/** Import a creation script, then split it back into one script per table. */
async function splitPerTable(script: string, varNames: string[]): Promise<PerTableResult> {
  const e = makeEditor();
  try {
    const g = buildCreationScriptGraph(script);
    await applyGraphToEditor(g, e.flow);
    return emitCreationScriptsForTables(e.flow, varNames);
  } finally {
    destroyEditor(e);
  }
}

const lines = (s: string): string[] => s.split('\n').filter((l) => l.trim() !== '');

const TIMESTAMP_RE = /\s*\/\/\{"timestamp":\s*(\d+)\}\s*$/;
/** Lines with the trailing `//{"timestamp": …}` comment stripped. */
const stripTs = (s: string): string[] => lines(s).map((l) => l.replace(TIMESTAMP_RE, ''));
/** Parsed timestamps, one per non-empty line. */
const timestamps = (s: string): number[] =>
  lines(s).map((l) => Number(TIMESTAMP_RE.exec(l)?.[1] ?? NaN));

category('Flow: creation script emit', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('producer → assignment, no leaked optional args', async () => {
    const r = await roundTrip('T = OpenFile("System:AppData/x.csv")');
    expect(r.warnings.length, 0, r.warnings.join(' ; '));
    // Optional sheetName seeded by the importer must NOT leak.
    expect(r.script, 'T = OpenFile("System:AppData/x.csv")');
  });

  test('in-place mutators become bare calls in execution order', async () => {
    const r = await roundTrip([
      'T = OpenFile("p.csv")',
      'AddNewColumn(T, "1", "a")',
      'AddNewColumn(T, "2", "b")',
    ].join('\n'));
    expect(r.warnings.length, 0, r.warnings.join(' ; '));
    const ls = lines(r.script);
    expect(ls[0], 'T = OpenFile("p.csv")', 'producer first, as an assignment');
    expect(ls[1], 'AddNewColumn(T, "1", "a")', 'first mutator is a bare call');
    expect(ls[2], 'AddNewColumn(T, "2", "b")', 'second mutator after the first');
  });

  test('variable references are bare identifiers, never quoted', async () => {
    const r = await roundTrip('T = OpenFile("p.csv")\nAddNewColumn(T, "1", "a")');
    expect(r.script.includes('AddNewColumn(T,'), true, 'table arg is the bare variable T');
    expect(r.script.includes('AddNewColumn("T"'), false, 'never the quoted name');
  });

  test('join: tables-by-name, column lists, default joinType omitted, inPlace kept', async () => {
    const r = await roundTrip(
      'Result = JoinTables("demog", "demog (2)", ["USUBJID"], ["USUBJID"], ' +
      '["USUBJID", "AGE"], ["USUBJID", "AGE"], "inner", true)');
    expect(r.warnings.length, 0, r.warnings.join(' ; '));
    expect(r.script.startsWith('Result = JoinTables("demog", "demog (2)"'), true, 'tables as quoted names');
    expect(r.script.includes('["USUBJID", "AGE"]'), true, 'column lists as arrays');
    expect(r.script.includes('inPlace = true'), true, 'non-default optional kept as a named arg');
    expect(r.script.includes('"inner"'), false, 'default joinType omitted');
  });

  test('friendly-name table reference resolves; order edge adds no line', async () => {
    const r = await roundTrip([
      'Mol1KLocal = OpenFile("local.csv")',
      'Result = JoinTables("mol1K local", "demog", ["prID"], ["prID"], ["prID"], ["prID"])',
    ].join('\n'));
    const ls = lines(r.script);
    expect(ls.length, 2, `expected exactly 2 lines, got: ${r.script}`);
    expect(ls[0], 'Mol1KLocal = OpenFile("local.csv")');
    expect(ls[1].startsWith('Result = JoinTables("mol1K local", "demog"'), true, 'table referenced by name');
  });

  test('emit → import → emit is stable (idempotent)', async () => {
    const script = [
      'T = OpenFile("p.csv")',
      'AddNewColumn(T, "1", "a")',
      'AddNewColumn(T, "2", "b")',
    ].join('\n');
    const first = (await roundTrip(script)).script;
    const second = (await roundTrip(first)).script;
    expect(second, first, 'a second round-trip reproduces the first');
  });

  test('an Output node anchors the producer exactly like a SetVar', async () => {
    const openFunc = DG.Func.find({name: 'OpenFile'})[0];
    if (!openFunc) {
      expect(true, true); // no live backend — nothing to check
      return;
    }
    const e = makeEditor();
    try {
      const open = await addNode(e.flow, ensureFuncNodeType(openFunc));
      open.inputValues['fullPath'] = 'p.csv';
      const out = await addNode(e.flow, 'Outputs/Table Output', 240, 0);
      out.properties['paramName'] = 'T';
      const outKey = Object.keys(open.outputs).find((k) => !isExecKey(k) && !k.endsWith('__pt'))!;
      await e.flow.addConnectionByKeys(open.id, outKey, out.id, 'table');

      const r = emitCreationScript(e.flow);
      expect(r.warnings.length, 0, r.warnings.join(' ; '));
      expect(r.script, 'T = OpenFile("p.csv")',
        'one anchored assignment — no intermediate variable, no redundant T = T line');
    } finally {
      destroyEditor(e);
    }
  });

  test('unsupported JS-only node warns and is skipped, rest still emits', async () => {
    const e = makeEditor();
    try {
      const s = await addNode(e.flow, 'Constants/String', 0, 0);
      s.properties['value'] = 'hello';
      const c = await addNode(e.flow, 'Comparisons/Contains', 240, 0);
      await e.flow.addConnectionByKeys(s.id, 'value', c.id, 'text');

      const r = emitCreationScript(e.flow);
      expect(r.warnings.some((w) => w.includes('no creation-script equivalent')), true,
        `expected an unsupported-node warning, got: ${r.warnings.join(' ; ')}`);
    } finally {
      destroyEditor(e);
    }
  });

  test('full chem script round-trips faithfully', async () => {
    if (!chemAvailable()) {
      expect(true, true);
      return;
    }
    const r = await roundTrip([
      'Mol1KLocal = OpenFile("System:AppData/Chem/mol1K.csv")',
      'Mol1K = OpenFile("System:AppData/Chem/mol1K.csv")',
      'Chem:addChemPropertiesColumns(Mol1K, "molecule", true, true, true, true, true, false, false, false, false)',
      'AddNewColumn(Mol1K, "${LogP} + ${MW} + ${HBD}", "smth", subscribeOnChanges = true)',
      'Result = JoinTables("mol1K", "mol1K local", ["prID"], ["prID"], ' +
        '["molecule", "prID", "smth"], ["molecule", "prID"])',
      'Chem:addChemPropertiesColumns(Result, "molecule", false, false, false, false, ' +
        'false, false, false, true, false)',
      'Chem:addChemPropertiesColumns(Result, "mol1K local.molecule", false, false, ' +
        'false, false, false, false, true, false, false)',
    ].join('\n'));
    expect(r.warnings.length, 0, r.warnings.join(' ; '));

    const s = r.script;
    expect(s.includes('Mol1K = OpenFile("System:AppData/Chem/mol1K.csv")'), true, 'producer assignment');
    expect(s.includes('Chem:addChemPropertiesColumns(Mol1K, "molecule", true,'), true, 'namespaced bare mutator');
    // GROK-20428 flipped subscribeOnChanges' default to true, so the platform
    // serializer omits it (optionals equal to their default are dropped) —
    // assert the call itself, agnostic of whether the named arg is present.
    expect(s.includes('AddNewColumn(Mol1K, "${LogP} + ${MW} + ${HBD}", "smth"'),
      true, 'AddNewColumn re-emits as a bare mutator call');
    expect(s.includes('Result = JoinTables("mol1K", "mol1K local"'), true, 'join produces Result');
    expect(s.includes('"mol1K local.molecule"'), true, 'qualified column name preserved');
    expect(s.includes('ResolveColumn'), false, 'no ResolveColumn leaks');
    expect(s.includes('.col('), false, 'columns are names, not JS table.col(...)');
  });

  // ---------- per-table split (Save Creation Scripts) ----------

  const TWO_TABLES = [
    'A = OpenFile("a.csv")',
    'AddNewColumn(A, "1", "x")',
    'B = OpenFile("b.csv")',
    'AddNewColumn(B, "2", "y")',
  ].join('\n');

  test('per-table: splits a combined script by owning variable', async () => {
    const r = await splitPerTable(TWO_TABLES, ['A', 'B']);
    expect(r.warnings.length, 0, r.warnings.join(' ; '));
    expect(r.unassigned.length, 0, `unexpected unassigned lines: ${r.unassigned.join(' ; ')}`);
    expect(r.tables.length, 2);

    expect(r.tables[0].variableName, 'A');
    expect(stripTs(r.tables[0].script).join('\n'),
      ['A = OpenFile("a.csv")', 'AddNewColumn(A, "1", "x")'].join('\n'),
      'table A keeps only its producer + mutator');

    expect(r.tables[1].variableName, 'B');
    expect(stripTs(r.tables[1].script).join('\n'),
      ['B = OpenFile("b.csv")', 'AddNewColumn(B, "2", "y")'].join('\n'),
      'table B keeps only its producer + mutator');
  });

  test('per-table: result is aligned to the requested order', async () => {
    const r = await splitPerTable(TWO_TABLES, ['B', 'A']);
    expect(r.tables.map((t) => t.variableName).join(','), 'B,A', 'order follows the request');
    expect(stripTs(r.tables[0].script)[0], 'B = OpenFile("b.csv")', 'first tab is B');
  });

  test('per-table: every line carries a strictly increasing timestamp comment', async () => {
    const r = await splitPerTable(TWO_TABLES, ['A', 'B']);
    for (const t of r.tables) {
      const ls = lines(t.script);
      expect(ls.length, 2, `${t.variableName}: two lines`);
      for (const l of ls)
        expect(TIMESTAMP_RE.test(l), true, `line lacks a timestamp comment: ${l}`);
      const ts = timestamps(t.script);
      expect(ts.every((n) => Number.isFinite(n)), true, 'all timestamps parse');
      for (let i = 1; i < ts.length; i++)
        expect(ts[i] > ts[i - 1], true, `timestamps must increase: ${ts.join(', ')}`);
    }
  });

  test('per-table: lines for an unrequested table fall into unassigned', async () => {
    const r = await splitPerTable(TWO_TABLES, ['A']);
    expect(r.tables.length, 1);
    expect(r.tables[0].variableName, 'A');
    expect(stripTs(r.tables[0].script).length, 2, 'A still fully built');
    // B's producer + mutator belong to no requested table.
    expect(r.unassigned.length, 2, `expected B's 2 lines unassigned, got: ${r.unassigned.join(' ; ')}`);
    expect(r.unassigned.some((l) => l.includes('B = OpenFile("b.csv")')), true);
  });

  test('per-table: a requested table absent from the graph yields an empty script', async () => {
    const r = await splitPerTable(TWO_TABLES, ['A', 'B', 'Missing']);
    expect(r.tables.length, 3);
    expect(r.tables[2].variableName, 'Missing');
    expect(r.tables[2].script, '', 'no lines own the missing variable');
  });
});
