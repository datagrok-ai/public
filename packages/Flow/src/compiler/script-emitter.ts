import {LGraph} from 'litegraph.js';
import {CompiledStep, compileGraph} from './graph-compiler';
import {FuncFlowNode} from '../types/funcflow-node';

export interface ScriptSettings {
  name: string;
  description: string;
  tags: string[];
}

export interface EmitOptions {
  /** When true, wraps each step in try/catch with event-firing instrumentation */
  instrumented?: boolean;
  /** UUID for this run — required when instrumented=true */
  runId?: string;
  /** Enable breakpoint node code emission (debug mode) */
  enableBreakpoints?: boolean;
  /** Halt execution on first error (default true) */
  haltOnError?: boolean;
}

/** Generates a valid Datagrok script from compiled steps.
 * When options.instrumented is true, wraps each step in try/catch with
 * event-firing code for live execution visualization. */
export function emitScript(graph: LGraph, settings: ScriptSettings, options?: EmitOptions): string {
  const steps = compileGraph(graph);
  const lines: string[] = [];
  const inst = options?.instrumented === true;

  // --- Header ---
  lines.push(`//name: ${settings.name}`);
  if (settings.description)
    lines.push(`//description: ${settings.description}`);
  lines.push('//language: javascript');
  if (settings.tags.length > 0)
    lines.push(`//tags: ${settings.tags.join(', ')}`);

  // --- Input declarations (table/dataframe inputs first) ---
  const inputSteps = steps.filter((s) => s.nodeType === 'input');
  inputSteps.sort((a, b) => {
    const aNode = graph.getNodeById(a.nodeId) as FuncFlowNode | undefined;
    const bNode = graph.getNodeById(b.nodeId) as FuncFlowNode | undefined;
    const aIsTable = aNode?.dgOutputType === 'dataframe';
    const bIsTable = bNode?.dgOutputType === 'dataframe';
    if (aIsTable && !bIsTable) return -1;
    if (!aIsTable && bIsTable) return 1;
    return 0;
  });
  for (const step of inputSteps) {
    const node = graph.getNodeById(step.nodeId) as FuncFlowNode | undefined;
    if (!node) continue;
    const inputLine = buildInputLine(step, node);
    if (inputLine) lines.push(inputLine);
  }

  // --- Output declarations ---
  const outputSteps = steps.filter((s) => s.nodeType === 'output');
  for (const step of outputSteps) {
    const node = graph.getNodeById(step.nodeId) as FuncFlowNode | undefined;
    if (!node) continue;
    const outputLine = buildOutputLine(step, node);
    if (outputLine) lines.push(outputLine);
  }

  lines.push('');

  // --- Instrumentation preamble ---
  if (inst)
    lines.push(...emitPreamble(options!.runId!));

  // --- Body ---
  for (const step of steps) {
    if (step.nodeType === 'input') continue;

    if (step.nodeType === 'output') {
      const inputExpr = step.inputs.get('value') || 'undefined';
      if (inputExpr !== 'undefined') {
        if (inst) {
          lines.push(...wrapInstrumented(
            `${step.variableName} = ${inputExpr};`,
            step, options!, {outputExpr: step.variableName},
          ));
        } else
          lines.push(`${step.variableName} = ${inputExpr};`);
      }
      continue;
    }

    if (step.nodeType === 'utility') {
      // Handle breakpoint nodes
      if (step.funcName === 'Breakpoint') {
        if (inst && options?.enableBreakpoints)
          lines.push(...emitBreakpointCode(step));
        continue;
      }
      const codeLine = emitUtilityStep(step);
      if (codeLine) {
        if (inst)
          lines.push(...wrapInstrumented(codeLine, step, options!, {outputExpr: step.variableName}));
        else
          lines.push(codeLine);
      }
      continue;
    }

    // Function step
    if (inst)
      lines.push(...emitFuncStepInstrumented(step, options!));
    else
      lines.push(emitFuncStep(step));
  }

  if (inst)
    lines.push(`__ff_emit('run-complete', -1, {success: true});`);

  return lines.join('\n');
}

function buildInputLine(step: CompiledStep, node: FuncFlowNode): string | null {
  const paramName = step.properties['paramName'];
  const outputType = node.dgOutputType || 'dynamic';
  const description = step.properties['description'] || '';

  let line = `//input: ${outputType} ${paramName}`;

  const qualifiers: string[] = [];

  const defaultVal = step.properties['defaultValue'];
  if (defaultVal !== undefined && defaultVal !== '' && defaultVal !== null)
    line = `//input: ${outputType} ${paramName} = ${formatHeaderDefault(defaultVal, outputType)}`;

  if (step.properties['typeFilter'])
    qualifiers.push(`type: ${step.properties['typeFilter']}`);
  if (step.properties['semTypeFilter'])
    qualifiers.push(`semType: ${step.properties['semTypeFilter']}`);
  if (step.properties['semType'])
    qualifiers.push(`semType: ${step.properties['semType']}`);
  if (step.properties['nullable'] === true)
    qualifiers.push('nullable: true');
  if (step.properties['caption'])
    qualifiers.push(`caption: ${step.properties['caption']}`);
  if (step.properties['choices']) {
    const items = step.properties['choices'].split(',').map((s: string) => s.trim()).filter(Boolean);
    qualifiers.push(`choices: [${items.map((s: string) => `"${s}"`).join(', ')}]`);
  }
  if (step.properties['min'] !== undefined && step.properties['min'] !== '')
    qualifiers.push(`min: ${step.properties['min']}`);
  if (step.properties['max'] !== undefined && step.properties['max'] !== '')
    qualifiers.push(`max: ${step.properties['max']}`);
  if (step.properties['showSlider'] === true)
    qualifiers.push('showSlider: true');

  if (qualifiers.length > 0)
    line += ` {${qualifiers.join('; ')}}`;

  if (description)
    line += ` [${description}]`;

  return line;
}

function buildOutputLine(step: CompiledStep, node: FuncFlowNode): string | null {
  const paramName = step.properties['paramName'];
  // For Value Output nodes, use the user-selected outputType property;
  // for Table Output, use the fixed dgOutputType ('dataframe')
  const outputType = step.properties['outputType'] || node.dgOutputType || 'dynamic';
  return `//output: ${outputType} ${paramName}`;
}

function emitFuncStep(step: CompiledStep): string {
  const params: string[] = [];
  for (const [name, expr] of step.inputs.entries())
    params.push(`${name}: ${expr}`);
  const paramsStr = params.length > 0 ? `{${params.join(', ')}}` : '{}';
  return `let ${step.variableName} = await grok.functions.call('${step.funcName}', ${paramsStr});`;
}

function emitUtilityStep(step: CompiledStep): string | null {
  switch (step.funcName) {
  case 'Select Column': {
    const tableExpr = step.inputs.get('table') || 'undefined';
    const colName = step.properties['columnName'] || '';
    return `let ${step.variableName} = ${tableExpr}.col('${colName}');`;
  }
  case 'Select Columns': {
    const tableExpr = step.inputs.get('table') || 'undefined';
    const colNames = (step.properties['columnNames'] || '').split(',').map((n: string) => n.trim()).filter(Boolean);
    const colExprs = colNames.map((n: string) => `${tableExpr}.col('${n}')`).join(', ');
    return `let ${step.variableName} = [${colExprs}];`;
  }
  case 'Add Table View': {
    const tableExpr = step.inputs.get('table') || 'undefined';
    return `let ${step.variableName} = grok.shell.addTableView(${tableExpr});`;
  }
  case 'Log': {
    const valueExpr = step.inputs.get('value') || 'undefined';
    const label = step.properties['label'] || '';
    return label ?
      `console.log('${label}:', ${valueExpr});` :
      `console.log(${valueExpr});`;
  }
  case 'Info': {
    const msgExpr = step.inputs.get('message') || '\'\'';
    return `grok.shell.info(${msgExpr});`;
  }
  case 'Warning': {
    const msgExpr = step.inputs.get('message') || '\'\'';
    return `grok.shell.warning(${msgExpr});`;
  }
  case 'ToString': {
    const valueExpr = step.inputs.get('value') || 'undefined';
    return `let ${step.variableName} = (${valueExpr}).toString();`;
  }
  case 'String': {
    const val = step.properties['value'] ?? '';
    return `let ${step.variableName} = ${JSON.stringify(String(val))};`;
  }
  case 'Int': {
    const val = Math.round(Number(step.properties['value'] ?? 0));
    return `let ${step.variableName} = ${val};`;
  }
  case 'Double': {
    const val = Number(step.properties['value'] ?? 0);
    return `let ${step.variableName} = ${val};`;
  }
  case 'Boolean': {
    const val = step.properties['value'] ? 'true' : 'false';
    return `let ${step.variableName} = ${val};`;
  }
  case 'List': {
    const raw = String(step.properties['value'] ?? '');
    const items = raw.split(',').map((s: string) => s.trim()).filter(Boolean);
    const arrStr = items.map((s: string) => JSON.stringify(s)).join(', ');
    return `let ${step.variableName} = [${arrStr}];`;
  }
  case 'Equals (==)': {
    const a = step.inputs.get('a') || 'undefined';
    const b = step.inputs.get('b') || 'undefined';
    return `let ${step.variableName} = (${a}) == (${b});`;
  }
  case 'Not Equals (!=)': {
    const a = step.inputs.get('a') || 'undefined';
    const b = step.inputs.get('b') || 'undefined';
    return `let ${step.variableName} = (${a}) != (${b});`;
  }
  case 'Greater Than (>)': {
    const a = step.inputs.get('a') || 'undefined';
    const b = step.inputs.get('b') || 'undefined';
    return `let ${step.variableName} = (${a}) > (${b});`;
  }
  case 'Greater Or Equal (>=)': {
    const a = step.inputs.get('a') || 'undefined';
    const b = step.inputs.get('b') || 'undefined';
    return `let ${step.variableName} = (${a}) >= (${b});`;
  }
  case 'Less Than (<)': {
    const a = step.inputs.get('a') || 'undefined';
    const b = step.inputs.get('b') || 'undefined';
    return `let ${step.variableName} = (${a}) < (${b});`;
  }
  case 'Less Or Equal (<=)': {
    const a = step.inputs.get('a') || 'undefined';
    const b = step.inputs.get('b') || 'undefined';
    return `let ${step.variableName} = (${a}) <= (${b});`;
  }
  case 'Contains': {
    const text = step.inputs.get('text') || '\'\'';
    const substring = step.inputs.get('substring') || '\'\'';
    return `let ${step.variableName} = (${text}).includes(${substring});`;
  }
  case 'Starts With': {
    const text = step.inputs.get('text') || '\'\'';
    const prefix = step.inputs.get('prefix') || '\'\'';
    return `let ${step.variableName} = (${text}).startsWith(${prefix});`;
  }
  case 'Ends With': {
    const text = step.inputs.get('text') || '\'\'';
    const suffix = step.inputs.get('suffix') || '\'\'';
    return `let ${step.variableName} = (${text}).endsWith(${suffix});`;
  }
  case 'Is Null': {
    const val = step.inputs.get('value') || 'undefined';
    return `let ${step.variableName} = (${val}) == null;`;
  }
  case 'FromJSON': {
    const jsonExpr = step.inputs.get('json') || '\'{}\'';
    return `let ${step.variableName} = JSON.parse(${jsonExpr});`;
  }
  case 'ToJSON': {
    const valExpr = step.inputs.get('value') || 'undefined';
    return `let ${step.variableName} = JSON.stringify(${valExpr});`;
  }
  default:
    return `// Unknown utility: ${step.funcName}`;
  }
}

// --- Instrumentation helpers ---

function emitPreamble(runId: string): string[] {
  return [
    `const __ff_runId = '${runId}';`,
    `const __ff_ch = 'funcflow.exec.' + __ff_runId;`,
    'function __ff_summarize(v) {',
    '  if (v == null) return {type:\'null\', value:null};',
    '  if (v.rowCount !== undefined && v.columns !== undefined)',
    '    return {type:\'dataframe\', rows:v.rowCount, cols:v.columns.length,',
    '      colNames:v.columns.names(), clone:v.clone()};',
    '  if (v.length !== undefined && v.name !== undefined && v.toList)',
    '    return {type:\'column\', name:v.name, length:v.length, sample:v.toList().slice(0,5)};',
    '  if (typeof v === \'object\') return {type:\'object\', str:String(v).slice(0,200)};',
    '  return {type:\'primitive\', value:v};',
    '}',
    'function __ff_emit(type, nodeId, data) {',
    '  grok.events.fireCustomEvent(__ff_ch, Object.assign({type, nodeId, timestamp:Date.now()}, data||{}));',
    '}',
    `__ff_emit('run-start', -1);`,
    '',
  ];
}

function wrapInstrumented(
  codeLine: string, step: CompiledStep, options: EmitOptions,
  extra?: {outputExpr?: string},
): string[] {
  const lines: string[] = [];
  lines.push(`__ff_emit('node-start', ${step.nodeId});`);

  // Hoist variable declaration before try block so it's accessible to downstream nodes
  let bodyLine = codeLine;
  const letMatch = codeLine.match(/^let\s+(\w+)\s*=/);
  if (letMatch) {
    lines.push(`let ${letMatch[1]};`);
    bodyLine = codeLine.replace(/^let\s+/, '');
  }

  lines.push('try {');
  lines.push(`  ${bodyLine}`);
  if (extra?.outputExpr) {
    lines.push(`  __ff_emit('node-complete', ${step.nodeId}, ` +
      `{outputs:{${extra.outputExpr}: __ff_summarize(${extra.outputExpr})}});`);
  } else
    lines.push(`  __ff_emit('node-complete', ${step.nodeId});`);
  lines.push('} catch (__ff_err) {');
  lines.push(`  __ff_emit('node-error', ${step.nodeId}, {error: __ff_err.message, stack: __ff_err.stack});`);
  if (options.haltOnError !== false) {
    lines.push(`  __ff_emit('run-complete', -1, {success: false});`);
    lines.push('  throw __ff_err;');
  }
  lines.push('}');
  return lines;
}

function emitFuncStepInstrumented(step: CompiledStep, options: EmitOptions): string[] {
  const params: string[] = [];
  for (const [name, expr] of step.inputs.entries())
    params.push(`${name}: ${expr}`);
  const paramsStr = params.length > 0 ? `{${params.join(', ')}}` : '{}';

  const lines: string[] = [];
  lines.push(`__ff_emit('node-start', ${step.nodeId});`);
  lines.push(`let ${step.variableName};`);
  lines.push('try {');
  lines.push(`  ${step.variableName} = await grok.functions.call('${step.funcName}', ${paramsStr});`);

  // Build outputs summary
  const outputEntries: string[] = [];
  for (const [, varName] of step.outputs.entries())
    outputEntries.push(`${varName}: __ff_summarize(${varName})`);
  const outputsObj = outputEntries.length > 0 ? `{${outputEntries.join(', ')}}` : '{}';
  lines.push(`  __ff_emit('node-complete', ${step.nodeId}, {outputs: ${outputsObj}});`);

  lines.push('} catch (__ff_err) {');
  lines.push(`  __ff_emit('node-error', ${step.nodeId}, {error: __ff_err.message, stack: __ff_err.stack});`);
  if (options.haltOnError !== false) {
    lines.push(`  __ff_emit('run-complete', -1, {success: false});`);
    lines.push('  throw __ff_err;');
  }
  lines.push('}');
  return lines;
}

function emitBreakpointCode(step: CompiledStep): string[] {
  return [
    `__ff_emit('breakpoint-hit', ${step.nodeId});`,
    'await new Promise((resolve) => {',
    `  const sub = grok.events.onCustomEvent(__ff_ch + '.continue').subscribe(() => {`,
    '    sub.unsubscribe();',
    '    resolve();',
    '  });',
    '});',
    `__ff_emit('node-complete', ${step.nodeId});`,
  ];
}

function formatHeaderDefault(value: any, type: string): string {
  switch (type) {
  case 'string':
    return `"${String(value)}"`;
  case 'bool':
    return value ? 'true' : 'false';
  case 'int':
  case 'double':
  case 'num':
    return String(value);
  default:
    return String(value);
  }
}
