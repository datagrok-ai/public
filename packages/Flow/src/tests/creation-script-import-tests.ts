import * as DG from 'datagrok-api/dg';
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes, registerAllFunctions} from '../rete/node-factory';
import {
  buildCreationScriptGraph, applyGraphToEditor, BuiltGraph,
  estimateNodeWidth, estimateNodeHeight,
} from '../import/creation-script-importer';
import {emitScript} from '../compiler/script-emitter';
import {validateGraph} from '../compiler/validator';
import {
  makeEditor, destroyEditor, nodesByFunc, oneNodeByFunc, nodesByLabel, sourceOf, until,
} from './test-utils';

const PASSTHROUGH = '__pt';
const SETTINGS = {name: 'Imported', description: '', tags: ['funcflow']};

/** The exact example from the feature request (uses Chem). */
const CHEM_PROPS_CALL =
  'Chem:addChemPropertiesColumns(Mol1K, "molecule", true, true, true, true, false, false, false, false, false)';
const CHEM_SCRIPT = [
  'Mol1K = OpenFile("System:AppData/Chem/mol1K.csv") //{"timestamp": 1781179268730}',
  `${CHEM_PROPS_CALL} //{"timestamp": 1781179278878}`,
  'AddNewColumn(Mol1K, "${HBA}+${HBD}+${LogP}", "sumOfSome", subscribeOnChanges = true) //{"timestamp": 1781179299560}',
].join('\n');

function chemAvailable(): boolean {
  try {
    return DG.Func.find({name: 'addChemPropertiesColumns'}).length > 0;
  } catch {
    return false;
  }
}

/** The SetVar node registering the given variable (the terminal per variable). */
function setVarFor(graph: BuiltGraph, varName: string) {
  return nodesByFunc(graph, 'SetVar').find((n) => n.inputValues['variableName'] === varName);
}

category('Flow: creation script import', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('simple assignment builds func node + SetVar, no output node', async () => {
    const g = buildCreationScriptGraph('T = OpenFile("System:AppData/x.csv")');
    expect(g.outputVariables.length, 1);
    expect(g.outputVariables[0], 'T');
    expect(g.nodes.some((n) => n.dgNodeType === 'output'), false, 'no output nodes');
    const open = oneNodeByFunc(g, 'OpenFile');
    expect(open.inputValues['fullPath'], 'System:AppData/x.csv'); // string slot → inputValue
    const setVar = setVarFor(g, 'T');
    expect(setVar != null, true, 'SetVar created');
    expect(setVar!.label, 'set: T');
    expect(sourceOf(g, setVar!, 'value')?.node, open, 'SetVar wired to OpenFile');
  });

  test('every variable gets its own SetVar node', async () => {
    const g = buildCreationScriptGraph([
      'A = OpenFile("a.csv")',
      'B = OpenFile("b.csv")',
      'AddNewColumn(B, "1", "x")',
    ].join('\n'));
    expect(g.outputVariables.join(','), 'A,B');
    const setVars = nodesByFunc(g, 'SetVar');
    expect(setVars.length, 2);
    const setA = setVarFor(g, 'A');
    const setB = setVarFor(g, 'B');
    expect(setA != null && setB != null, true, 'one SetVar per variable');

    // A is untouched → wired straight to its OpenFile result.
    const opens = nodesByFunc(g, 'OpenFile');
    const openA = opens.find((n) => n.inputValues['fullPath'] === 'a.csv')!;
    expect(sourceOf(g, setA!, 'value')?.node, openA);
    // B was mutated → wired to AddNewColumn's pass-through (final state).
    const add = oneNodeByFunc(g, 'AddNewColumn');
    const srcB = sourceOf(g, setB!, 'value');
    expect(srcB?.node, add);
    expect(srcB!.key.endsWith(PASSTHROUGH), true);
  });

  test('imported nodes are collapsed by default', async () => {
    const g = buildCreationScriptGraph('T = OpenFile("p.csv")\nAddNewColumn(T, "1", "a")');
    expect(g.nodes.length > 0, true);
    for (const node of g.nodes)
      expect(node.collapsed, true, `node "${node.label}" should be collapsed`);
  });

  test('primitive on editable slot → inputValue; non-editable handled separately', async () => {
    const g = buildCreationScriptGraph('T = OpenFile("p.csv")\nAddNewColumn(T, "1+2", "c", subscribeOnChanges = true)');
    const add = oneNodeByFunc(g, 'AddNewColumn');
    expect(add.inputValues['expression'], '1+2'); // string slot
    expect(add.inputValues['name'], 'c'); // string slot
    expect(add.inputValues['subscribeOnChanges'], true); // bool slot
  });

  test('bare mutating calls thread the table through pass-through outputs (execution order)', async () => {
    const g = buildCreationScriptGraph([
      'T = OpenFile("p.csv")',
      'AddNewColumn(T, "1", "a")',
      'AddNewColumn(T, "2", "b")',
    ].join('\n'));
    const adds = nodesByFunc(g, 'AddNewColumn');
    expect(adds.length, 2);
    const [first, second] = adds;
    // second AddNewColumn's table comes from the first's pass-through, not OpenFile.
    const tableSrc = sourceOf(g, second, 'table');
    expect(tableSrc?.node, first, 'second consumes first');
    expect(tableSrc!.key.endsWith(PASSTHROUGH), true, 'via pass-through');
    // SetVar stores the final state of the table → second AddNewColumn's pass-through.
    const setVar = setVarFor(g, 'T')!;
    const setSrc = sourceOf(g, setVar, 'value');
    expect(setSrc?.node, second);
    expect(setSrc!.key.endsWith(PASSTHROUGH), true);
  });

  test('trailing // metadata comments and blank lines are stripped', async () => {
    const g = buildCreationScriptGraph([
      'T = OpenFile("System:AppData/x.csv") //{"timestamp": 123}',
      '',
      'AddNewColumn(T, "1", "a") //{"timestamp": 456}',
    ].join('\n'));
    expect(g.warnings.length, 0, `unexpected warnings: ${g.warnings.join(' | ')}`);
    expect(nodesByFunc(g, 'OpenFile').length, 1);
    expect(nodesByFunc(g, 'AddNewColumn').length, 1);
  });

  test('URL with // inside string survives comment stripping', async () => {
    const g = buildCreationScriptGraph('T = OpenFile("https://example.com/data.csv") //{"timestamp": 1}');
    const open = oneNodeByFunc(g, 'OpenFile');
    expect(open.inputValues['fullPath'], 'https://example.com/data.csv');
  });

  test('empty script throws', async () => {
    let threw = false;
    try {
      buildCreationScriptGraph('\n// only comments\n');
    } catch {
      threw = true;
    }
    expect(threw, true);
  });

  // ---------- the exact reported bug (Chem) ----------

  test('column argument becomes an inline column input value, no Select Column node', async () => {
    if (!chemAvailable()) return; // Chem not on this server — skip gracefully
    const g = buildCreationScriptGraph(CHEM_SCRIPT);
    expect(g.warnings.length, 0, `unexpected warnings: ${g.warnings.join(' | ')}`);

    const open = oneNodeByFunc(g, 'OpenFile');
    const chem = oneNodeByFunc(g, 'addChemPropertiesColumns');

    // ResolveColumn (broken on the platform) is no longer wired through a
    // Select Column node — the column name is stored on the input instead.
    expect(nodesByFunc(g, 'ResolveColumn').length, 0, 'no ResolveColumn node');
    expect(nodesByLabel(g, 'Select Column').length, 0, 'no Select Column node');
    expect(chem.inputValues['molecules'], 'molecule', 'column name stored as input value');
    expect(sourceOf(g, chem, 'molecules'), null, 'molecules input is not connected');

    // The column resolves against the chem call's own table input.
    expect((chem.properties['columnTables'] as Record<string, string>)['molecules'], 'table');

    // The chem call's own inputs.
    expect(sourceOf(g, chem, 'table')?.node, open);
    expect(chem.inputValues['MW'], true); // boolean slot → inputValue
    expect(chem.inputValues['logS'], false);
  });

  test('full chem script: ordering through chem + AddNewColumn to output', async () => {
    if (!chemAvailable()) return;
    const g = buildCreationScriptGraph(CHEM_SCRIPT);
    const chem = oneNodeByFunc(g, 'addChemPropertiesColumns');
    const add = oneNodeByFunc(g, 'AddNewColumn');

    // AddNewColumn runs after chem: its table comes from chem's pass-through.
    const addTableSrc = sourceOf(g, add, 'table');
    expect(addTableSrc?.node, chem, 'AddNewColumn ordered after chem');
    expect(addTableSrc!.key.endsWith(PASSTHROUGH), true);

    // SetVar(Mol1K) stores the final Mol1K = AddNewColumn pass-through.
    const setVar = setVarFor(g, 'Mol1K')!;
    const setSrc = sourceOf(g, setVar, 'value');
    expect(setSrc?.node, add);
    expect(g.outputVariables.join(','), 'Mol1K');
  });

  // ---------- integration: apply to a live editor + emit ----------

  test('applied graph validates and emits ordered script', async () => {
    const e = makeEditor();
    try {
      const g = buildCreationScriptGraph([
        'T = OpenFile("System:AppData/x.csv")',
        'AddNewColumn(T, "1", "a")',
        'AddNewColumn(T, "2", "b")',
      ].join('\n'));
      await applyGraphToEditor(g, e.flow);

      expect(e.flow.getNodeCount(), g.nodes.length);
      expect(e.flow.getConnectionCount(), g.connections.length);

      const errors = validateGraph(e.flow).filter((r) => r.severity === 'error');
      expect(errors.length, 0, `validation errors: ${errors.map((x) => x.message).join('; ')}`);

      const script = emitScript(e.flow, SETTINGS);
      const iOpen = script.indexOf('OpenFile');
      const iAdd = script.indexOf('AddNewColumn');
      expect(iOpen >= 0 && iAdd >= 0, true, 'both calls present');
      expect(iOpen < iAdd, true, 'OpenFile emitted before AddNewColumn');
    } finally {
      destroyEditor(e);
    }
  });

  test('column_list arguments become inline comma-separated input values', async () => {
    // Each column_list parses to an array of ResolveColumn calls; the names are
    // joined into the input value, and the column→table association pairs by the
    // numeric suffix (keys2/values2 → table2).
    const g = buildCreationScriptGraph(
      'Result = JoinTables("demog", "demog (2)", ["USUBJID"], ["USUBJID"], ' +
      '["USUBJID", "AGE", "SEX"], ["USUBJID", "AGE"], "inner", true)');
    expect(g.warnings.length, 0, `unexpected warnings: ${g.warnings.join(' | ')}`);

    const join = oneNodeByFunc(g, 'JoinTables');
    expect(join.inputValues['joinType'], 'inner');
    expect(join.inputValues['inPlace'], true);

    // Table name strings parse to ResolveTable calls — substituted with the
    // Select Table utility (grok.shell.tableByName), titled after the table.
    expect(nodesByFunc(g, 'ResolveTable').length, 0, 'no ResolveTable nodes');
    const table1 = sourceOf(g, join, 'table1')!.node;
    const table2 = sourceOf(g, join, 'table2')!.node;
    expect(table1.dgTypeName, 'Utilities/Select Table');
    expect(table2.dgTypeName, 'Utilities/Select Table');
    expect(table1.properties['tableName'], 'demog');
    expect(table2.properties['tableName'], 'demog (2)');

    // Column lists are inlined as comma-separated input values — no Select Columns nodes.
    expect(nodesByLabel(g, 'Select Columns').length, 0, 'no Select Columns nodes');
    expect(join.inputValues['keys1'], 'USUBJID');
    expect(join.inputValues['keys2'], 'USUBJID');
    expect(join.inputValues['values1'], 'USUBJID, AGE, SEX');
    expect(join.inputValues['values2'], 'USUBJID, AGE');

    // Numbered pairing recorded in the association: *1 lists → table1, *2 → table2.
    const assoc = join.properties['columnTables'] as Record<string, string>;
    expect(assoc['keys1'], 'table1');
    expect(assoc['values1'], 'table1');
    expect(assoc['keys2'], 'table2');
    expect(assoc['values2'], 'table2');

    // The variable Result is stored by SetVar from JoinTables' real output.
    const setVar = setVarFor(g, 'Result')!;
    const setSrc = sourceOf(g, setVar, 'value');
    expect(setSrc?.node, join);
    expect(setSrc!.key.endsWith(PASSTHROUGH), false, 'stored from the real result, not a pass-through');

    // Exact graph size: JoinTables + 2 Select Table + SetVar (no output, no Select Columns).
    expect(g.nodes.length, 4);
  });

  test('layout: a producer path sits above the path that consumes its table', async () => {
    // The Join (defined FIRST) reads "Second" via a Select Table; the Second
    // producer path is defined LATER. Layout must still place the producer band
    // above the consumer band (dependency order beats script order), so the
    // top-first execution order resolves the table before the join uses it.
    const g = buildCreationScriptGraph([
      'Joined = JoinTables("First", "Second", ["Id"], ["Id"], ["Id"], ["Id"])',
      'Second = OpenFile("s.csv")',
      'AddNewColumn(Second, "1", "x")',
    ].join('\n'));

    const join = oneNodeByFunc(g, 'JoinTables');
    const setSecond = setVarFor(g, 'Second')!;
    const openSecond = nodesByFunc(g, 'OpenFile')[0];

    // Bands are vertically disjoint, so every Second-path node is above the join.
    expect(openSecond.pos.y < join.pos.y, true, 'producer OpenFile above the join');
    expect(setSecond.pos.y < join.pos.y, true, 'producer SetVar above the join');

    // Sanity: the Select Table the join reads is named "Second".
    const tables = nodesByFunc(g, 'JoinTables').length; // ensure exactly one join
    expect(tables, 1);
    const t2 = sourceOf(g, join, 'table2')!.node;
    expect(t2.properties['tableName'], 'Second');
  });

  test('layout: edges point right, no node boxes overlap', async () => {
    const g = buildCreationScriptGraph([
      'Result = JoinTables("demog", "demog (2)", ["USUBJID"], ["USUBJID"], ' +
        '["USUBJID", "AGE", "SEX"], ["USUBJID", "AGE"], "inner", true)',
      'AddNewColumn(Result, "1", "x")',
      'Other = OpenFile("o.csv")',
    ].join('\n'));

    // Direction: every connection flows left → right.
    for (const c of g.connections) {
      expect(c.source.pos.x < c.target.pos.x, true,
        `"${c.source.label}" (${c.source.pos.x}) must be left of "${c.target.label}" (${c.target.pos.x})`);
    }

    // No overlap between estimated bounding boxes.
    for (let i = 0; i < g.nodes.length; i++) {
      for (let j = i + 1; j < g.nodes.length; j++) {
        const a = g.nodes[i];
        const b = g.nodes[j];
        const separated =
          a.pos.x + estimateNodeWidth(a) <= b.pos.x || b.pos.x + estimateNodeWidth(b) <= a.pos.x ||
          a.pos.y + estimateNodeHeight(a) <= b.pos.y || b.pos.y + estimateNodeHeight(b) <= a.pos.y;
        expect(separated, true, `"${a.label}" and "${b.label}" overlap`);
      }
    }
  });

  test('connections render when nodes start collapsed', async () => {
    const e = makeEditor();
    try {
      const g = buildCreationScriptGraph('T = OpenFile("p.csv")\nAddNewColumn(T, "1", "a")');
      await applyGraphToEditor(g, e.flow);

      // Collapsed nodes render socket DOM only for connected sockets; the
      // editor must re-render endpoints when connections arrive after the
      // node (regression: wires were invisible until expand + collapse).
      const rendered = await until(() => {
        const sockets = e.container.querySelectorAll('.ff-node-collapsed-sockets .ff-socket');
        const paths = Array.from(e.container.querySelectorAll('.ff-connection-path'));
        return sockets.length > 0 && paths.length === g.connections.length &&
          paths.every((p) => (p.getAttribute('d') ?? '') !== '');
      });
      expect(rendered, true, 'collapsed sockets and connection paths must render without expanding');
    } finally {
      destroyEditor(e);
    }
  });

  test('chem script emits table.col(...) instead of ResolveColumn', async () => {
    if (!chemAvailable()) return;
    const e = makeEditor();
    try {
      const g = buildCreationScriptGraph(CHEM_SCRIPT);
      await applyGraphToEditor(g, e.flow);

      const errors = validateGraph(e.flow).filter((r) => r.severity === 'error');
      expect(errors.length, 0, `validation errors: ${errors.map((x) => x.message).join('; ')}`);

      const script = emitScript(e.flow, SETTINGS);
      expect(script.includes(`.col('molecule')`), true, 'column selected via table.col()');
      expect(script.includes('ResolveColumn'), false, 'no ResolveColumn in generated script');
      // Execution order preserved: OpenFile → addChemPropertiesColumns → AddNewColumn.
      const iOpen = script.indexOf('OpenFile');
      const iChem = script.indexOf('addChemPropertiesColumns');
      const iAdd = script.indexOf('AddNewColumn');
      expect(iOpen < iChem && iChem < iAdd, true, 'calls emitted in script order');
    } finally {
      destroyEditor(e);
    }
  });

  test('SetVar surfaces its value for the output preview (instrumented run)', async () => {
    const e = makeEditor();
    try {
      const g = buildCreationScriptGraph('T = OpenFile("System:AppData/x.csv")');
      await applyGraphToEditor(g, e.flow);
      const script = emitScript(e.flow, SETTINGS, {instrumented: true, runId: 'run-1'});
      // SetVar's node-complete event captures the value under the variable name,
      // so clicking the node opens the docked preview (table → grid, etc.) even
      // though SetVar declares no output.
      expect(script.includes('"T": __ff_summarize('), true, 'SetVar value surfaced for preview');
    } finally {
      destroyEditor(e);
    }
  });

  test('SetVar emits a runtime-guarded registration under the dataframe name', async () => {
    const e = makeEditor();
    try {
      const g = buildCreationScriptGraph('T = OpenFile("System:AppData/x.csv")');
      await applyGraphToEditor(g, e.flow);

      // The value slot can be `dynamic`, so the dataframe check is done at
      // runtime: when the value is a DataFrame, SetVar also registers it under
      // the table's runtime .name, so GetVars that use the actual table name
      // resolve (single node on the canvas, two assignments in the output).
      const script = emitScript(e.flow, SETTINGS);
      expect(script.includes('instanceof DG.DataFrame'), true, 'runtime dataframe guard emitted');
      expect(/variableName: \w+\.name\b/.test(script), true, 'second SetVar keyed by the dataframe runtime name');
      // The variable-name registration is still present (single node, two assigns).
      expect(script.includes(`variableName: "T"`), true, 'primary SetVar by variable name');
    } finally {
      destroyEditor(e);
    }
  });

  test('mixed local/connected script: columns inlined, no Select Column(s) nodes, compiles', async () => {
    if (!chemAvailable()) return;
    const e = makeEditor();
    try {
      const FULL_SCRIPT = [
        'Mol1KLocal = OpenTable("65d4d9d0-48b0-11f1-e424-4b91b3dfc6ce")',
        'Mol1K = OpenFile("System:AppData/Chem/mol1K.csv") // {"timestamp": 1781796656926}',
        'Chem:addChemPropertiesColumns(Mol1K, "molecule", true, true, true, true, true, false, false, false, false)',
        'AddNewColumn(Mol1K, "${LogP} + ${MW} + ${HBD}", "smth", subscribeOnChanges = true)',
        'Result = JoinTables("mol1K", "mol1K local", ["prID"], ["prID"], ' +
          '["molecule", "prID", "IDDB", "MW", "HBA", "HBD", "LogP", "LogS", "smth"], ' +
          '["molecule", "prID", "IDDB"])',
        'Chem:addChemPropertiesColumns(Result, "molecule", false, false, false, false, ' +
          'false, false, false, true, false)',
        'Chem:addChemPropertiesColumns(Result, "mol1K local.molecule", false, false, ' +
          'false, false, false, false, true, false, false)',
      ].join('\n');

      const g = buildCreationScriptGraph(FULL_SCRIPT);
      await applyGraphToEditor(g, e.flow);

      // The whole point: no Select Column / Select Columns clutter.
      expect(nodesByLabel(g, 'Select Column').length, 0, 'no Select Column nodes');
      expect(nodesByLabel(g, 'Select Columns').length, 0, 'no Select Columns nodes');

      const errors = validateGraph(e.flow).filter((r) => r.severity === 'error');
      expect(errors.length, 0, `validation errors: ${errors.map((x) => x.message).join('; ')}`);

      const script = emitScript(e.flow, SETTINGS);
      expect(script.includes('ResolveColumn'), false, 'no ResolveColumn in generated script');
      // Single-column args inlined as table.col(...), incl. the qualified name.
      expect(script.includes(`.col('molecule')`), true, 'molecule column selected via .col()');
      expect(script.includes(`.col('mol1K local.molecule')`), true, 'qualified column name preserved');
      // column_list args inlined as arrays of table.col(...).
      expect(script.includes(`.col('prID')`), true, 'join key column via .col()');
      expect(script.includes(`.col('smth')`), true, 'join value column via .col()');
    } finally {
      destroyEditor(e);
    }
  });

  test('join script emits grok.shell.tableByName instead of ResolveTable', async () => {
    const e = makeEditor();
    try {
      const g = buildCreationScriptGraph(
        'Result = JoinTables("demog", "demog (2)", ["USUBJID"], ["USUBJID"], ' +
        '["USUBJID", "AGE"], ["USUBJID", "AGE"], "inner", true)');
      await applyGraphToEditor(g, e.flow);

      const errors = validateGraph(e.flow).filter((r) => r.severity === 'error');
      expect(errors.length, 0, `validation errors: ${errors.map((x) => x.message).join('; ')}`);

      const script = emitScript(e.flow, SETTINGS);
      expect(script.includes('grok.shell.tableByName("demog")'), true, 'table1 via tableByName');
      expect(script.includes('grok.shell.tableByName("demog (2)")'), true, 'table2 via tableByName');
      expect(script.includes('ResolveTable'), false, 'no ResolveTable in generated script');
      expect(script.includes(`.col('USUBJID')`), true, 'key columns via table.col()');
    } finally {
      destroyEditor(e);
    }
  });
});
