import * as DG from 'datagrok-api/dg';
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes, registerAllFunctions} from '../rete/node-factory';
import {buildCreationScriptGraph, applyGraphToEditor} from '../import/creation-script-importer';
import {emitCreationScript, CreationScriptResult} from '../compiler/creation-script-emitter';
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

const lines = (s: string): string[] => s.split('\n').filter((l) => l.trim() !== '');

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
    expect(s.includes('AddNewColumn(Mol1K, "${LogP} + ${MW} + ${HBD}", "smth", subscribeOnChanges = true)'),
      true, 'optional bool as a named arg');
    expect(s.includes('Result = JoinTables("mol1K", "mol1K local"'), true, 'join produces Result');
    expect(s.includes('"mol1K local.molecule"'), true, 'qualified column name preserved');
    expect(s.includes('ResolveColumn'), false, 'no ResolveColumn leaks');
    expect(s.includes('.col('), false, 'columns are names, not JS table.col(...)');
  });
});
