import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes, registerAllFunctions, ensureFuncNodeType, getRegisteredFuncs} from '../rete/node-factory';
import {topologicalSort} from '../compiler/topological-sort';
import {emitScript} from '../compiler/script-emitter';
import {emitCreationScript} from '../compiler/creation-script-emitter';
import {validateGraph} from '../compiler/validator';
import {ExecutionController} from '../execution/execution-controller';
import {NodeExecStatus} from '../execution/execution-state';
import {makeEditor, destroyEditor, addNode, until} from './test-utils';

const SETTINGS = {name: 'TestFlow', description: 'test', tags: ['funcflow']};

category('Flow: topological sort', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('orders source before target', async () => {
    const e = makeEditor();
    try {
      const input = await addNode(e.flow, 'Inputs/Table Input');
      const output = await addNode(e.flow, 'Outputs/Table Output');
      await e.flow.addConnectionByKeys(input.id, 'table', output.id, 'table');
      const sorted = topologicalSort(e.flow);
      expect(sorted.length, 2);
      expect(sorted.indexOf(input.id) < sorted.indexOf(output.id), true);
    } finally {
      destroyEditor(e);
    }
  });

  test('chained pass-through preserves order', async () => {
    const e = makeEditor();
    try {
      const c = await addNode(e.flow, 'Constants/String');
      const t1 = await addNode(e.flow, 'Utilities/ToString');
      const t2 = await addNode(e.flow, 'Utilities/ToString');
      await e.flow.addConnectionByKeys(c.id, 'value', t1.id, 'value');
      await e.flow.addConnectionByKeys(t1.id, 'text', t2.id, 'value');
      const sorted = topologicalSort(e.flow);
      expect(sorted.indexOf(c.id) < sorted.indexOf(t1.id), true);
      expect(sorted.indexOf(t1.id) < sorted.indexOf(t2.id), true);
    } finally {
      destroyEditor(e);
    }
  });

  test('disjoint subgraphs execute top-first, fully', async () => {
    const e = makeEditor();
    try {
      // Bottom chain created FIRST — canvas position must win over creation order.
      const bottomIn = await addNode(e.flow, 'Constants/String', 0, 500);
      const bottomOut = await addNode(e.flow, 'Utilities/ToString', 300, 500);
      await e.flow.addConnectionByKeys(bottomIn.id, 'value', bottomOut.id, 'value');
      const topIn = await addNode(e.flow, 'Constants/String', 0, 50);
      const topOut = await addNode(e.flow, 'Utilities/ToString', 300, 50);
      await e.flow.addConnectionByKeys(topIn.id, 'value', topOut.id, 'value');

      const sorted = topologicalSort(e.flow);
      const at = (id: string): number => sorted.indexOf(id);
      expect(Math.max(at(topIn.id), at(topOut.id)) < Math.min(at(bottomIn.id), at(bottomOut.id)), true,
        'top chain must fully precede the bottom chain');
    } finally {
      destroyEditor(e);
    }
  });

  test('within a component, ready nodes process top-to-bottom', async () => {
    const e = makeEditor();
    try {
      // the lower source is created first — order must come from position, not creation
      const lower = await addNode(e.flow, 'Constants/String', 0, 400);
      const upper = await addNode(e.flow, 'Constants/String', 0, 100);
      const sink = await addNode(e.flow, 'Comparisons/Equals (==)', 300, 250);
      await e.flow.addConnectionByKeys(upper.id, 'value', sink.id, 'a');
      await e.flow.addConnectionByKeys(lower.id, 'value', sink.id, 'b');

      const sorted = topologicalSort(e.flow);
      expect(sorted.indexOf(upper.id) < sorted.indexOf(lower.id), true, 'upper source first');
      expect(sorted.indexOf(sink.id), sorted.length - 1, 'sink last');
    } finally {
      destroyEditor(e);
    }
  });

  test('detects cycles', async () => {
    const e = makeEditor();
    try {
      const a = await addNode(e.flow, 'Utilities/ToString');
      const b = await addNode(e.flow, 'Utilities/ToString');
      await e.flow.addConnectionByKeys(a.id, 'text', b.id, 'value');
      await e.flow.addConnectionByKeys(b.id, 'text', a.id, 'value');
      let threw = false;
      try {
        topologicalSort(e.flow);
      } catch {
        threw = true;
      }
      expect(threw, true);
    } finally {
      destroyEditor(e);
    }
  });
});

category('Flow: script emitter', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('emits input/output headers and pass-through body', async () => {
    const e = makeEditor();
    try {
      const input = await addNode(e.flow, 'Inputs/Table Input');
      input.properties['paramName'] = 'myTable';
      const output = await addNode(e.flow, 'Outputs/Table Output');
      output.properties['paramName'] = 'res';
      await e.flow.addConnectionByKeys(input.id, 'table', output.id, 'table');

      const script = emitScript(e.flow, SETTINGS);
      expect(script.includes('//input: dataframe myTable'), true, 'input header');
      expect(script.includes('//output: dataframe res'), true, 'output header');
      expect(script.includes('res = myTable;'), true, 'output assignment');
      expect(script.includes('//language: javascript'), true, 'language header');
    } finally {
      destroyEditor(e);
    }
  });

  test('constant string is inlined into Value Output', async () => {
    const e = makeEditor();
    try {
      const c = await addNode(e.flow, 'Constants/String');
      c.properties['value'] = 'hello';
      // emission must dispatch on the registered type, not the user-editable label
      c.label = 'const: hello';
      const out = await addNode(e.flow, 'Outputs/Value Output');
      out.properties['paramName'] = 'greeting';
      await e.flow.addConnectionByKeys(c.id, 'value', out.id, 'value');

      const script = emitScript(e.flow, SETTINGS);
      expect(script.includes('"hello"'), true, 'string literal present');
      expect(/greeting\s*=\s*\w+;/.test(script), true, 'output assigned from upstream variable');
    } finally {
      destroyEditor(e);
    }
  });

  test('instrumented mode wraps steps with run events', async () => {
    const e = makeEditor();
    try {
      const input = await addNode(e.flow, 'Inputs/Table Input');
      const output = await addNode(e.flow, 'Outputs/Table Output');
      await e.flow.addConnectionByKeys(input.id, 'table', output.id, 'table');

      const script = emitScript(e.flow, SETTINGS, {instrumented: true, runId: 'run-123'});
      expect(script.includes('run-123'), true, 'run id embedded');
      expect(script.includes('run-complete'), true, 'completion event emitted');
    } finally {
      destroyEditor(e);
    }
  });

  test('OpenFile: sheetName stays out of the call unless the path is xlsx', async () => {
    const info = getRegisteredFuncs().find((f) => f.func.name === 'OpenFile');
    if (!info) return;
    const e = makeEditor();
    try {
      const node = await addNode(e.flow, info.nodeTypeName);
      node.inputValues['fullPath'] = 'System:AppData/Chem/mol1K.sdf';
      let script = emitScript(e.flow, SETTINGS);
      // OpenFile forwards any non-null sheetName as a second importer argument — the sdf importer takes one
      expect(/grok\.functions\.call\('OpenFile', \{fullPath: "[^"]+"\}\)/.test(script), true,
        `sdf call carries only the path (script: ${script})`);
      expect(script.includes('sheetName'), false, 'no sheetName for a non-xlsx path');

      node.inputValues['sheetName'] = 'Sheet1';
      script = emitScript(e.flow, SETTINGS);
      expect(script.includes('sheetName'), false, 'a typed sheet name is still dropped for sdf');

      node.inputValues['fullPath'] = 'System:AppData/Demo/book.xlsx';
      script = emitScript(e.flow, SETTINGS);
      expect(script.includes('sheetName: "Sheet1"'), true, `xlsx keeps the sheet name (script: ${script})`);
    } finally {
      destroyEditor(e);
    }
  });

  test('an untouched optional scalar is omitted so the function default applies', async () => {
    const info = getRegisteredFuncs().find((f) => f.func.name === 'OpenFile');
    if (!info) return;
    const e = makeEditor();
    try {
      const node = await addNode(e.flow, info.nodeTypeName);
      node.inputValues['fullPath'] = 'a.csv';
      const script = emitScript(e.flow, SETTINGS);
      expect(script.includes('sheetName'), false, 'blank optional omitted');
      expect(script.includes('fullPath: "a.csv"'), true, 'filled value kept');
    } finally {
      destroyEditor(e);
    }
  });

  test('an Open File node pointed at an sdf runs end-to-end', async () => {
    const info = getRegisteredFuncs().find((f) => f.func.name === 'OpenFile');
    if (!info) return;
    const e = makeEditor();
    try {
      const node = await addNode(e.flow, info.nodeTypeName);
      node.inputValues['fullPath'] = 'System:AppData/Chem/mol1K.sdf';
      const ctrl = new ExecutionController(e.flow);
      expect(ctrl.runAutorun(new Set(), SETTINGS), 'started', 'run starts');
      const done = (): boolean => {
        const s = ctrl.state.getNodeState(node.id)?.status;
        return s === NodeExecStatus.completed || s === NodeExecStatus.errored;
      };
      expect(await until(done, 30000), true, 'run finished');
      const st = ctrl.state.getNodeState(node.id);
      expect(st?.status, NodeExecStatus.completed, `sdf opened via the emitted call (${st?.error ?? ''})`);
    } finally {
      destroyEditor(e);
    }
  });

  test('side-effect-only utilities (Log) emit no output summary when instrumented', async () => {
    const e = makeEditor();
    try {
      const c = await addNode(e.flow, 'Constants/String');
      c.properties['value'] = 'hello';
      const log = await addNode(e.flow, 'Utilities/Log');
      await e.flow.addConnectionByKeys(c.id, 'value', log.id, 'value');
      const toStr = await addNode(e.flow, 'Utilities/ToString');
      await e.flow.addConnectionByKeys(c.id, 'value', toStr.id, 'value');

      const script = emitScript(e.flow, SETTINGS, {instrumented: true, runId: 'run-log'});
      expect(script.includes('console.log('), true, 'the Log body is emitted');
      expect(script.includes('__ff_summarize(log'), false, 'no summary of an undeclared variable');
      expect(script.includes(`__ff_stash(${JSON.stringify(log.id)}`), false,
        'no stash of an undeclared variable (the compiler gives Log no phantom output)');
      expect(script.includes(`__ff_emit('node-complete', '${log.id}');`), true,
        'the Log node completes with no outputs payload');
      // instrumented `let` is hoisted out of the try block, so the body line is a bare assignment
      const toStringVar = script.match(/(\w+) = \(.*\)\.toString\(\);/)?.[1];
      expect(!!toStringVar, true, 'the ToString step declares its variable');
      expect(script.includes(`__ff_summarize(${toStringVar})`), true, 'ToString keeps its output summary');
    } finally {
      destroyEditor(e);
    }
  });

  test('Select Table fails fast when no table resolves', async () => {
    const e = makeEditor();
    try {
      const sel = await addNode(e.flow, 'Utilities/Select Table');
      sel.properties['tableName'] = 'My Table';
      const out = await addNode(e.flow, 'Outputs/Table Output', 240, 0);
      out.properties['paramName'] = 'res';
      await e.flow.addConnectionByKeys(sel.id, 'table', out.id, 'table');

      const script = emitScript(e.flow, SETTINGS);
      expect(script.includes('grok.shell.tableByName("My Table")'), true, 'name lookup emitted');
      expect(script.includes('throw new Error("Select Table: no open table or variable named'), true,
        'a null table throws with the table name instead of failing downstream');
      const inst = emitScript(e.flow, SETTINGS, {instrumented: true, runId: 'run-sel'});
      expect(inst.includes('throw new Error("Select Table: no open table or variable named'), true,
        'the guard survives instrumentation');
    } finally {
      destroyEditor(e);
    }
  });

  test('Output and SetVar compile to the same contract', async () => {
    const setVarFunc = DG.Func.find({name: 'SetVar'})[0];
    if (!setVarFunc) {
      expect(true, true); // no live backend — nothing to check
      return;
    }
    const e = makeEditor();
    try {
      const input = await addNode(e.flow, 'Inputs/Table Input');
      input.properties['paramName'] = 'myTable';
      const output = await addNode(e.flow, 'Outputs/Table Output', 240, 0);
      output.properties['paramName'] = 'res';
      await e.flow.addConnectionByKeys(input.id, 'table', output.id, 'table');
      const setVar = await addNode(e.flow, ensureFuncNodeType(setVarFunc), 240, 200);
      setVar.inputValues['variableName'] = 'MyResult';
      await e.flow.addConnectionByKeys(input.id, 'table', setVar.id, 'value');

      const script = emitScript(e.flow, SETTINGS);
      expect(script.includes(
        `await grok.functions.call('SetVar', {variableName: "res", value: myTable});`), true,
      'an Output node doubles as a SetVar registration');
      expect(script.includes('//output: dataframe MyResult'), true,
        'a SetVar node declares an output header with the connection-inferred type');
      expect(script.includes('MyResult = myTable;'), true,
        'a SetVar node assigns its output variable');
    } finally {
      destroyEditor(e);
    }
  });

  test('a no-output mutator with two table inputs captures both modified tables', async () => {
    const script = DG.Script.create([
      '//name: TwoTableMutator',
      '//language: javascript',
      '//input: dataframe t1',
      '//input: dataframe t2',
      'grok.shell.info("noop");',
    ].join('\n'));
    expect(script.outputs.length, 0, 'the synthetic mutator declares no outputs');

    const typeName = ensureFuncNodeType(script);
    const e = makeEditor();
    try {
      const in1 = await addNode(e.flow, 'Inputs/Table Input', 0, 0);
      const in2 = await addNode(e.flow, 'Inputs/Table Input', 0, 140);
      in1.properties['paramName'] = 'a';
      in2.properties['paramName'] = 'b';
      const mut = await addNode(e.flow, typeName, 320, 60);
      await e.flow.addConnectionByKeys(in1.id, 'table', mut.id, 't1');
      await e.flow.addConnectionByKeys(in2.id, 'table', mut.id, 't2');

      const inst = emitScript(e.flow, SETTINGS, {instrumented: true, runId: 'r1'});
      expect(inst.includes(`'t1 (modified)': __ff_summarize(`), true, 'first input table captured');
      expect(inst.includes(`'t2 (modified)': __ff_summarize(`), true, 'second input table captured');
    } finally {
      destroyEditor(e);
    }
  });

  test('dataframe pass-throughs get dims-only summaries; utility summaries are slot-keyed', async () => {
    const script = DG.Script.create([
      '//name: TableThreader',
      '//language: javascript',
      '//input: dataframe t',
      '//output: int result',
      'result = t.rowCount;',
    ].join('\n'));
    const typeName = ensureFuncNodeType(script);
    const e = makeEditor();
    try {
      const src = await addNode(e.flow, 'Inputs/Table Input', 0, 0);
      const fn = await addNode(e.flow, typeName, 320, 0);
      await e.flow.addConnectionByKeys(src.id, 'table', fn.id, 't');

      const inst = emitScript(e.flow, SETTINGS, {instrumented: true, runId: 'r1'});
      expect(inst.includes('function __ff_dims'), true, 'dims helper in the preamble');
      expect(inst.includes('"t__pt": __ff_dims('), true, 'pass-through summarized (dims-only)');

      const sel = await addNode(e.flow, 'Utilities/Select Table', 0, 200);
      sel.properties['tableName'] = 'demog';
      const inst2 = emitScript(e.flow, SETTINGS, {instrumented: true, runId: 'r2'});
      const selOutKey = Object.keys(sel.outputs).find((k) => !k.startsWith('__'))!;
      expect(inst2.includes(`{outputs:{${JSON.stringify(selOutKey)}: __ff_summarize(`), true,
        `utility summary keyed by slot key "${selOutKey}"`);
    } finally {
      destroyEditor(e);
    }
  });
});

category('Flow: multi-output funcs', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('reads outputs off the call object; variable names are valid identifiers', async () => {
    // the name starts with a digit and declares two outputs — both failure modes at once
    const script = DG.Script.create([
      '//name: 2InputsFlow',
      '//language: javascript',
      '//input: int firstNum',
      '//input: int secondNum',
      '//output: dataframe result1',
      '//output: dataframe result2',
      'result1 = grok.data.demo.demog();',
      'result2 = grok.data.demo.demog();',
    ].join('\n'));
    expect(script.outputs.length, 2, 'the synthetic func has two outputs');

    const typeName = ensureFuncNodeType(script);
    const e = makeEditor();
    try {
      const fn = await addNode(e.flow, typeName);
      const out1 = await addNode(e.flow, 'Outputs/Value Output', 340, 20);
      const out2 = await addNode(e.flow, 'Outputs/Value Output', 340, 140);
      out1.properties['paramName'] = 'Result';
      out2.properties['paramName'] = 'Result2';
      await e.flow.addConnectionByKeys(fn.id, 'result1', out1.id, 'value');
      await e.flow.addConnectionByKeys(fn.id, 'result2', out2.id, 'value');

      const clean = emitScript(e.flow, SETTINGS);

      expect(/\blet\s+[0-9]/.test(clean), false, 'no variable name starts with a digit');
      const m = clean.match(/let\s+(\w+)\s*=\s*await grok\.functions\.call\(/);
      expect(m != null, true, 'the func call is assigned to a variable');
      const callVar = m![1];
      expect(/^[0-9]/.test(callVar), false, `call var "${callVar}" must not start with a digit`);

      expect(clean.includes(`Result = ${callVar}.result1;`), true, 'Result assigned from the result1 property');
      expect(clean.includes(`Result2 = ${callVar}.result2;`), true, 'Result2 assigned from the result2 property');
      expect(clean.includes(`${callVar}_result1`), false, 'no reference to an undeclared per-output variable');
      expect(clean.includes(`if (${callVar}.result1 != null) await ${callVar}.result1.meta.detectSemanticTypes();`),
        true, 'semantic detection runs on the property expression');

      const inst = emitScript(e.flow, SETTINGS, {instrumented: true, runId: 'r1'});
      expect(inst.includes(`"result1": __ff_summarize(${callVar}.result1`), true,
        'summary keyed by slot key, value reads the property');
      expect(inst.includes(`"result2": __ff_summarize(${callVar}.result2`), true);
    } finally {
      destroyEditor(e);
    }
  });
});

category('Flow: validator', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('empty graph yields a warning', async () => {
    const e = makeEditor();
    try {
      const results = validateGraph(e.flow);
      expect(results.length >= 1, true);
      expect(results[0].severity, 'warning');
    } finally {
      destroyEditor(e);
    }
  });

  test('column input without table input is an error', async () => {
    const e = makeEditor();
    try {
      await addNode(e.flow, 'Inputs/Column Input');
      const results = validateGraph(e.flow);
      expect(results.some((r) => r.severity === 'error'), true);
    } finally {
      destroyEditor(e);
    }
  });

  test('valid input → output graph has no errors', async () => {
    const e = makeEditor();
    try {
      const input = await addNode(e.flow, 'Inputs/Table Input');
      const output = await addNode(e.flow, 'Outputs/Table Output');
      await e.flow.addConnectionByKeys(input.id, 'table', output.id, 'table');
      const results = validateGraph(e.flow);
      expect(results.some((r) => r.severity === 'error'), false);
    } finally {
      destroyEditor(e);
    }
  });

  test('duplicate variable names across outputs and SetVars are errors', async () => {
    const e = makeEditor();
    try {
      const input = await addNode(e.flow, 'Inputs/Table Input');
      input.properties['paramName'] = 'myTable';
      const out1 = await addNode(e.flow, 'Outputs/Table Output', 240, 0);
      out1.properties['paramName'] = 'res';
      const out2 = await addNode(e.flow, 'Outputs/Value Output', 240, 120);
      out2.properties['paramName'] = 'res';
      await e.flow.addConnectionByKeys(input.id, 'table', out1.id, 'table');
      await e.flow.addConnectionByKeys(input.id, 'table', out2.id, 'value');

      const isDup = (r: {severity: string; message: string}): boolean =>
        r.severity === 'error' && r.message.includes(`Duplicate variable name 'res'`);
      expect(validateGraph(e.flow).some(isDup), true, 'two outputs sharing a name is an error');

      out2.properties['paramName'] = 'other';
      expect(validateGraph(e.flow).some(isDup), false, 'distinct output names pass');

      const setVarFunc = DG.Func.find({name: 'SetVar'})[0];
      if (!setVarFunc) return; // no live backend — the output/output half is checked
      const sv = await addNode(e.flow, ensureFuncNodeType(setVarFunc), 240, 240);
      sv.inputValues['variableName'] = 'res';
      await e.flow.addConnectionByKeys(input.id, 'table', sv.id, 'value');
      expect(validateGraph(e.flow).some(isDup), true, 'a SetVar colliding with an output is an error');

      sv.inputValues['variableName'] = 'stored';
      expect(validateGraph(e.flow).some((r) => r.message.startsWith('Duplicate')), false,
        'distinct names across outputs and SetVars pass');
    } finally {
      destroyEditor(e);
    }
  });
});

category('Flow: func wrappers', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  const appendType = (): string | null =>
    getRegisteredFuncs().find((f) => f.func.name === 'AppendTables')?.nodeTypeName ?? null;

  test('a wrapped func node exposes the wrapper inputs, not the raw signature', async () => {
    const typeName = appendType();
    if (!typeName) return; // AppendTables not on this stand — skip
    const e = makeEditor();
    try {
      const node = await addNode(e.flow, typeName);
      expect(!!node.funcWrapper, true, 'wrapper read onto the node');
      expect('table1' in node.inputs, true, 'exposed table socket');
      expect('table2' in node.inputs, true, 'exposed table socket');
      expect('tables' in node.inputs, false, 'raw dataframe_list slot not exposed');
      expect('table1__pt' in node.outputs, true, 'pass-throughs mirror the exposed inputs');
      expect(node.passthroughCount, 2, 'pass-through count follows the wrapper');
      expect(node.requiredInputs.join(','), 'table1,table2', 'exposed tables gate the run');
      expect('result' in node.outputs, true, 'real output untouched');
    } finally {
      destroyEditor(e);
    }
  });

  test('compile folds the exposed inputs into the real call arguments', async () => {
    const typeName = appendType();
    if (!typeName) return;
    const e = makeEditor();
    try {
      const a = await addNode(e.flow, 'Inputs/Table Input', 0, 0);
      a.properties['paramName'] = 'tA';
      const b = await addNode(e.flow, 'Inputs/Table Input', 0, 200);
      b.properties['paramName'] = 'tB';
      const ap = await addNode(e.flow, typeName, 300, 100);
      await e.flow.addConnectionByKeys(a.id, 'table', ap.id, 'table1');
      await e.flow.addConnectionByKeys(b.id, 'table', ap.id, 'table2');

      const script = emitScript(e.flow, SETTINGS);
      expect(/grok\.functions\.call\('AppendTables', \{tables: \[\w+, \w+\]\}\)/.test(script), true,
        `real list argument emitted (script: ${script})`);
      expect(script.includes('table1:'), false, 'exposed names never leak into the call');

      const instrumented = emitScript(e.flow, SETTINGS, {instrumented: true, runId: 'w1'});
      expect(instrumented.includes('{tables: ['), true, 'instrumented call reshaped too');
      expect(instrumented.includes('"table1__pt":'), true,
        'pass-through stash still keyed by the exposed input');
    } finally {
      destroyEditor(e);
    }
  });

  test('the reshaped call shape is what the platform accepts', async () => {
    if (!appendType()) return;
    const dfA = DG.DataFrame.fromCsv('x\n1\n2');
    const dfB = DG.DataFrame.fromCsv('x\n3');
    const res = await grok.functions.call('AppendTables', {tables: [dfA, dfB]}) as DG.DataFrame;
    expect(res.rowCount, 3, 'a JS array of DataFrames marshals to the dataframe_list param');
  });

  test('a wrapped node has no creation-script equivalent — warns instead of emitting garbage', async () => {
    const typeName = appendType();
    if (!typeName) return;
    const e = makeEditor();
    try {
      await addNode(e.flow, typeName);
      const r = emitCreationScript(e.flow);
      expect(r.warnings.some((w) => w.includes('AppendTables')), true,
        `warning names the wrapped func (got: ${r.warnings.join(' ; ')})`);
      expect(r.script.includes('AppendTables'), false, 'no bogus call emitted');
    } finally {
      destroyEditor(e);
    }
  });
});
