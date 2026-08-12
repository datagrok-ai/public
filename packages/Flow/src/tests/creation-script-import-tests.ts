import * as DG from 'datagrok-api/dg';
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes, registerAllFunctions} from '../rete/node-factory';
import {EXEC_IN_KEY, EXEC_OUT_KEY} from '../rete/scheme';
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
    expect(open.inputValues['fullPath'], 'System:AppData/x.csv');
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

    const opens = nodesByFunc(g, 'OpenFile');
    const openA = opens.find((n) => n.inputValues['fullPath'] === 'a.csv')!;
    expect(sourceOf(g, setA!, 'value')?.node, openA);
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
    expect(add.inputValues['expression'], '1+2');
    expect(add.inputValues['name'], 'c');
    expect(add.inputValues['subscribeOnChanges'], true);
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
    const tableSrc = sourceOf(g, second, 'table');
    expect(tableSrc?.node, first, 'second consumes first');
    expect(tableSrc!.key.endsWith(PASSTHROUGH), true, 'via pass-through');
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

  test('column argument becomes an inline column input value, no Select Column node', async () => {
    if (!chemAvailable()) return; // Chem not on this server — skip gracefully
    const g = buildCreationScriptGraph(CHEM_SCRIPT);
    expect(g.warnings.length, 0, `unexpected warnings: ${g.warnings.join(' | ')}`);

    const open = oneNodeByFunc(g, 'OpenFile');
    const chem = oneNodeByFunc(g, 'addChemPropertiesColumns');

    expect(nodesByFunc(g, 'ResolveColumn').length, 0, 'no ResolveColumn node');
    expect(nodesByLabel(g, 'Select Column').length, 0, 'no Select Column node');
    expect(chem.inputValues['molecules'], 'molecule', 'column name stored as input value');
    expect(sourceOf(g, chem, 'molecules'), null, 'molecules input is not connected');

    expect((chem.properties['columnTables'] as Record<string, string>)['molecules'], 'table');

    expect(sourceOf(g, chem, 'table')?.node, open);
    expect(chem.inputValues['MW'], true);
    expect(chem.inputValues['logS'], false);
  });

  test('full chem script: ordering through chem + AddNewColumn to output', async () => {
    if (!chemAvailable()) return;
    const g = buildCreationScriptGraph(CHEM_SCRIPT);
    const chem = oneNodeByFunc(g, 'addChemPropertiesColumns');
    const add = oneNodeByFunc(g, 'AddNewColumn');

    const addTableSrc = sourceOf(g, add, 'table');
    expect(addTableSrc?.node, chem, 'AddNewColumn ordered after chem');
    expect(addTableSrc!.key.endsWith(PASSTHROUGH), true);

    const setVar = setVarFor(g, 'Mol1K')!;
    const setSrc = sourceOf(g, setVar, 'value');
    expect(setSrc?.node, add);
    expect(g.outputVariables.join(','), 'Mol1K');
  });

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
    const g = buildCreationScriptGraph(
      'Result = JoinTables("demog", "demog (2)", ["USUBJID"], ["USUBJID"], ' +
      '["USUBJID", "AGE", "SEX"], ["USUBJID", "AGE"], "inner", true)');
    expect(g.warnings.length, 0, `unexpected warnings: ${g.warnings.join(' | ')}`);

    const join = oneNodeByFunc(g, 'JoinTables');
    expect(join.inputValues['joinType'], 'inner');
    expect(join.inputValues['inPlace'], true);

    expect(nodesByFunc(g, 'ResolveTable').length, 0, 'no ResolveTable nodes');
    const table1 = sourceOf(g, join, 'table1')!.node;
    const table2 = sourceOf(g, join, 'table2')!.node;
    expect(table1.dgTypeName, 'Utilities/Select Table');
    expect(table2.dgTypeName, 'Utilities/Select Table');
    expect(table1.properties['tableName'], 'demog');
    expect(table2.properties['tableName'], 'demog (2)');

    expect(nodesByLabel(g, 'Select Columns').length, 0, 'no Select Columns nodes');
    expect(join.inputValues['keys1'], 'USUBJID');
    expect(join.inputValues['keys2'], 'USUBJID');
    expect(join.inputValues['values1'], 'USUBJID, AGE, SEX');
    expect(join.inputValues['values2'], 'USUBJID, AGE');

    const assoc = join.properties['columnTables'] as Record<string, string>;
    expect(assoc['keys1'], 'table1');
    expect(assoc['values1'], 'table1');
    expect(assoc['keys2'], 'table2');
    expect(assoc['values2'], 'table2');

    const setVar = setVarFor(g, 'Result')!;
    const setSrc = sourceOf(g, setVar, 'value');
    expect(setSrc?.node, join);
    expect(setSrc!.key.endsWith(PASSTHROUGH), false, 'stored from the real result, not a pass-through');

    // JoinTables + 2 Select Table + SetVar
    expect(g.nodes.length, 4);
  });

  test('layout: a producer path sits above the path that consumes its table', async () => {
    // the consumer join is defined FIRST — dependency order must beat script order
    const g = buildCreationScriptGraph([
      'Joined = JoinTables("First", "Second", ["Id"], ["Id"], ["Id"], ["Id"])',
      'Second = OpenFile("s.csv")',
      'AddNewColumn(Second, "1", "x")',
    ].join('\n'));

    const join = oneNodeByFunc(g, 'JoinTables');
    const setSecond = setVarFor(g, 'Second')!;
    const openSecond = nodesByFunc(g, 'OpenFile')[0];

    expect(openSecond.pos.y < join.pos.y, true, 'producer OpenFile above the join');
    expect(setSecond.pos.y < join.pos.y, true, 'producer SetVar above the join');

    const tables = nodesByFunc(g, 'JoinTables').length;
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

    for (const c of g.connections) {
      expect(c.source.pos.x < c.target.pos.x, true,
        `"${c.source.label}" (${c.source.pos.x}) must be left of "${c.target.label}" (${c.target.pos.x})`);
    }

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

      const script = emitScript(e.flow, SETTINGS);
      expect(script.includes('instanceof DG.DataFrame'), true, 'runtime dataframe guard emitted');
      expect(/variableName: \w+\.name\b/.test(script), true, 'second SetVar keyed by the dataframe runtime name');
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

      expect(nodesByLabel(g, 'Select Column').length, 0, 'no Select Column nodes');
      expect(nodesByLabel(g, 'Select Columns').length, 0, 'no Select Columns nodes');

      const errors = validateGraph(e.flow).filter((r) => r.severity === 'error');
      expect(errors.length, 0, `validation errors: ${errors.map((x) => x.message).join('; ')}`);

      const script = emitScript(e.flow, SETTINGS);
      expect(script.includes('ResolveColumn'), false, 'no ResolveColumn in generated script');
      expect(script.includes(`.col('molecule')`), true, 'molecule column selected via .col()');
      expect(script.includes(`.col('mol1K local.molecule')`), true, 'qualified column name preserved');
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

  function orderEdges(graph: BuiltGraph) {
    return graph.connections.filter((c) => c.order);
  }

  test('order edge inferred from a SetVar to a table referenced by friendly name', async () => {
    // "mol1K local" is Mol1KLocal's friendly name — no data edge to the producer,
    // so an order edge must force the producer first
    const g = buildCreationScriptGraph([
      'Mol1KLocal = OpenFile("local.csv")',
      'Result = JoinTables("mol1K local", "demog", ["prID"], ["prID"], ["prID"], ["prID"])',
    ].join('\n'));

    const edges = orderEdges(g);
    expect(edges.length, 1, 'exactly one order edge (only "mol1K local" matches a variable)');
    const edge = edges[0];

    const setLocal = setVarFor(g, 'Mol1KLocal')!;
    const selectLocal = g.nodes.find((n) =>
      n.dgTypeName === 'Utilities/Select Table' && n.properties['tableName'] === 'mol1K local')!;
    expect(edge.source, setLocal, 'order edge starts at the producing SetVar');
    expect(edge.target, selectLocal, 'order edge ends at the name-referenced Select Table');
    expect(edge.sourceKey, EXEC_OUT_KEY, 'wired from the exec-out port');
    expect(edge.targetKey, EXEC_IN_KEY, 'wired into the exec-in port');

    const selectDemog = g.nodes.find((n) =>
      n.dgTypeName === 'Utilities/Select Table' && n.properties['tableName'] === 'demog')!;
    expect(edges.some((c) => c.target === selectDemog), false, 'no order edge for the unmatched table');
  });

  test('no order edges when no table name matches a variable', async () => {
    const g = buildCreationScriptGraph([
      'A = OpenFile("a.csv")',
      'Result = JoinTables("totally different", "demog", ["id"], ["id"], ["id"], ["id"])',
    ].join('\n'));
    expect(orderEdges(g).length, 0, 'no inferred order edges');
  });

  test('inferred order edge applies to a live editor, sorts producer before the reference', async () => {
    const e = makeEditor();
    try {
      const g = buildCreationScriptGraph([
        'Mol1KLocal = OpenFile("local.csv")',
        'Result = JoinTables("mol1K local", "demog", ["prID"], ["prID"], ["prID"], ["prID"])',
      ].join('\n'));
      const added = await applyGraphToEditor(g, e.flow);
      expect(added, g.connections.length, 'every connection (incl. the order edge) applied');

      const errors = validateGraph(e.flow).filter((r) => r.severity === 'error');
      expect(errors.length, 0, `validation errors: ${errors.map((x) => x.message).join('; ')}`);

      const script = emitScript(e.flow, SETTINGS);
      const iSet = script.indexOf('"Mol1KLocal"');
      const iRef = script.indexOf('tableByName("mol1K local")');
      expect(iSet >= 0 && iRef >= 0, true, 'both the SetVar and the name reference are emitted');
      expect(iSet < iRef, true, 'Mol1KLocal is registered before it is read by name');
    } finally {
      destroyEditor(e);
    }
  });
});
