/** Generates the runnable Datagrok JavaScript from `CompiledStep[]`.
 *
 * Two emission modes:
 *   - clean: standard Datagrok script (with `//input:` / `//output:` headers)
 *   - instrumented: each step wrapped in try/catch firing `funcflow.exec.<runId>`
 *     events for live execution visualization. */

import {FlowEditor} from '../rete/flow-editor';
import {FlowNode} from '../rete/scheme';
import {CompiledStep, compileGraph} from './graph-compiler';

export interface ScriptSettings {
  name: string;
  description: string;
  tags: string[];
}

export interface EmitOptions {
  /** When true, wraps each step in try/catch with event-firing instrumentation. */
  instrumented?: boolean;
  /** UUID for this run — required when `instrumented=true`. */
  runId?: string;
  /** Whether to emit breakpoint pause code (debug mode). */
  enableBreakpoints?: boolean;
  /** Halt execution on first error (default true). */
  haltOnError?: boolean;
  /** When set, only compile/emit steps for these node ids (the rest are
   *  dropped). Used to run a slice up to a target node — the set must be
   *  closed under "ancestors" so no surviving step references a dropped one. */
  onlyNodeIds?: Set<string>;
}

const PASSTHROUGH_SUFFIX = '__pt';

export function emitScript(
  flow: FlowEditor, settings: ScriptSettings, options?: EmitOptions,
): string {
  let steps = compileGraph(flow);
  if (options?.onlyNodeIds) steps = steps.filter((s) => options.onlyNodeIds!.has(s.nodeId));
  const lines: string[] = [];
  const inst = options?.instrumented === true;

  // -------- header --------
  lines.push(`//name: ${settings.name}`);
  if (settings.description) lines.push(`//description: ${settings.description}`);
  lines.push('//language: javascript');
  if (settings.tags.length > 0) lines.push(`//tags: ${settings.tags.join(', ')}`);

  // -------- input declarations (dataframes first) --------
  const inputSteps = steps.filter((s) => s.nodeType === 'input');
  inputSteps.sort((a, b) => {
    const aIsTable = (flow.getNodeById(a.nodeId)?.dgOutputType) === 'dataframe';
    const bIsTable = (flow.getNodeById(b.nodeId)?.dgOutputType) === 'dataframe';
    return aIsTable === bIsTable ? 0 : aIsTable ? -1 : 1;
  });
  for (const step of inputSteps) {
    const node = flow.getNodeById(step.nodeId);
    if (!node) continue;
    const line = buildInputLine(step, node);
    if (line) lines.push(line);
  }

  // -------- output declarations --------
  for (const step of steps.filter((s) => s.nodeType === 'output')) {
    const node = flow.getNodeById(step.nodeId);
    if (!node) continue;
    lines.push(buildOutputLine(step, node));
  }

  lines.push('');

  if (inst) lines.push(...emitPreamble(options!.runId!));

  // -------- body --------
  for (const step of steps) {
    if (step.nodeType === 'input') continue;

    if (step.nodeType === 'output') {
      const inputKey = Array.from(step.inputs.keys())[0] ?? 'value';
      const inputExpr = step.inputs.get(inputKey) ?? 'undefined';
      if (inputExpr === 'undefined') continue;
      const node = flow.getNodeById(step.nodeId);
      const declaredType = (step.properties['outputType'] as string | undefined)
        ?? node?.dgOutputType ?? undefined;
      const line = `${step.variableName} = ${inputExpr};`;
      if (inst)
        lines.push(...wrapInstrumented(line, step, options!, {outputExpr: step.variableName, declaredType}));
      else
        lines.push(line);
      continue;
    }

    if (step.nodeType === 'utility') {
      if (step.funcName === 'Breakpoint') {
        if (inst && options?.enableBreakpoints) lines.push(...emitBreakpointCode(step));
        continue;
      }
      const code = emitUtilityStep(step);
      if (!code) continue;
      if (inst)
        lines.push(...wrapInstrumented(code, step, options!, {outputExpr: step.variableName}));
      else
        lines.push(code);
      continue;
    }

    // func step
    if (inst) lines.push(...emitFuncStepInstrumented(step, options!, flow));
    else lines.push(emitFuncStep(step));

    // A SetVar whose value is a dataframe at runtime also registers it under
    // the table's runtime name, so creation-script GetVars that reference the
    // table by its actual name resolve (the platform resolver registers both).
    // The dataframe check is emitted into the script (runtime instanceof), since
    // the value slot is often `dynamic` even when it carries a dataframe.
    const setVarValue = setVarValueExpr(step, flow);
    if (setVarValue) lines.push(emitSetVarByDataframeName(setVarValue));

    // Auto-detect semantic types on dataframe outputs.
    lines.push(...emitDetectSemanticTypes(step, flow));
  }

  if (inst) lines.push(`__ff_emit('run-complete', '', {success: true});`);

  return lines.join('\n');
}

function buildInputLine(step: CompiledStep, node: FlowNode): string | null {
  const paramName = String(step.properties['paramName']);
  const outputType = node.dgOutputType ?? 'dynamic';

  let line = `//input: ${outputType} ${paramName}`;
  const qualifiers: string[] = [];

  const defaultVal = step.properties['defaultValue'];
  if (defaultVal !== undefined && defaultVal !== '' && defaultVal !== null)
    line = `//input: ${outputType} ${paramName} = ${formatHeaderDefault(defaultVal, outputType)}`;

  if (step.properties['typeFilter']) qualifiers.push(`type: ${step.properties['typeFilter']}`);
  if (step.properties['semTypeFilter']) qualifiers.push(`semType: ${step.properties['semTypeFilter']}`);
  if (step.properties['semType']) qualifiers.push(`semType: ${step.properties['semType']}`);
  if (step.properties['nullable'] === true) qualifiers.push('nullable: true');
  if (step.properties['caption']) qualifiers.push(`caption: ${step.properties['caption']}`);
  if (step.properties['choices']) {
    const items = String(step.properties['choices']).split(',').map((s) => s.trim()).filter(Boolean);
    qualifiers.push(`choices: [${items.map((s) => `"${s}"`).join(', ')}]`);
  }
  if (step.properties['min'] !== undefined && step.properties['min'] !== '')
    qualifiers.push(`min: ${step.properties['min']}`);
  if (step.properties['max'] !== undefined && step.properties['max'] !== '')
    qualifiers.push(`max: ${step.properties['max']}`);
  if (step.properties['showSlider'] === true) qualifiers.push('showSlider: true');

  if (qualifiers.length > 0) line += ` {${qualifiers.join('; ')}}`;
  if (node.description) line += ` [${node.description}]`;
  return line;
}

function buildOutputLine(step: CompiledStep, node: FlowNode): string {
  const paramName = String(step.properties['paramName']);
  const outputType = (step.properties['outputType'] as string | undefined)
    ?? node.dgOutputType ?? 'dynamic';
  return `//output: ${outputType} ${paramName}`;
}

function emitFuncStep(step: CompiledStep): string {
  const params: string[] = [];
  for (const [name, expr] of step.inputs) params.push(`${name}: ${expr}`);
  const paramsStr = params.length > 0 ? `{${params.join(', ')}}` : '{}';
  return `let ${step.variableName} = await grok.functions.call('${step.funcName}', ${paramsStr});`;
}

// eslint-disable-next-line complexity
function emitUtilityStep(step: CompiledStep): string | null {
  switch (step.funcName) {
  case 'Select Column': {
    const t = step.inputs.get('table') ?? 'undefined';
    const n = step.properties['columnName'] ?? '';
    return `let ${step.variableName} = ${t}.col('${n}');`;
  }
  case 'Select Columns': {
    const t = step.inputs.get('table') ?? 'undefined';
    const names = String(step.properties['columnNames'] ?? '').split(',').map((s) => s.trim()).filter(Boolean);
    const exprs = names.map((n) => `${t}.col('${n}')`).join(', ');
    return `let ${step.variableName} = [${exprs}];`;
  }
  case 'Select Table': {
    // Resolve an open table by name, tolerating how it was registered: the
    // exact name, a no-spaces variant, and a lower-camel variant; for each,
    // try a shell table (tableByName) then a context variable (getVar). The
    // first non-null wins.
    const raw = String(step.properties['tableName'] ?? '');
    const noSpaces = raw.replace(/ /g, '');
    const lowerFirst = noSpaces.charAt(0).toLowerCase() + noSpaces.slice(1);
    const names = Array.from(new Set([raw, noSpaces, lowerFirst]));
    const expr = names.map((n) => {
      const q = JSON.stringify(n);
      return `grok.shell.tableByName(${q}) ?? grok.shell.getVar(${q})`;
    }).join(' ?? ');
    return `let ${step.variableName} = ${expr};`;
  }
  case 'Add Table View': {
    const t = step.inputs.get('table') ?? 'undefined';
    return `let ${step.variableName} = grok.shell.addTableView(${t});`;
  }
  case 'Log': {
    const v = step.inputs.get('value') ?? 'undefined';
    const label = step.properties['label'] ?? '';
    return label ? `console.log('${label}:', ${v});` : `console.log(${v});`;
  }
  case 'Info': {
    const m = step.inputs.get('message') ?? '\'\'';
    return `grok.shell.info(${m});`;
  }
  case 'Warning': {
    const m = step.inputs.get('message') ?? '\'\'';
    return `grok.shell.warning(${m});`;
  }
  case 'ToString': {
    const v = step.inputs.get('value') ?? 'undefined';
    return `let ${step.variableName} = (${v}).toString();`;
  }
  case 'String':
    return `let ${step.variableName} = ${JSON.stringify(String(step.properties['value'] ?? ''))};`;
  case 'Int':
    return `let ${step.variableName} = ${Math.round(Number(step.properties['value'] ?? 0))};`;
  case 'Double':
    return `let ${step.variableName} = ${Number(step.properties['value'] ?? 0)};`;
  case 'Boolean':
    return `let ${step.variableName} = ${step.properties['value'] ? 'true' : 'false'};`;
  case 'List': {
    const raw = String(step.properties['value'] ?? '');
    const items = raw.split(',').map((s) => s.trim()).filter(Boolean);
    return `let ${step.variableName} = [${items.map((s) => JSON.stringify(s)).join(', ')}];`;
  }
  case 'Equals (==)': {
    const a = step.inputs.get('a') ?? 'undefined'; const b = step.inputs.get('b') ?? 'undefined';
    return `let ${step.variableName} = (${a}) == (${b});`;
  }
  case 'Not Equals (!=)': {
    const a = step.inputs.get('a') ?? 'undefined'; const b = step.inputs.get('b') ?? 'undefined';
    return `let ${step.variableName} = (${a}) != (${b});`;
  }
  case 'Greater Than (>)': {
    const a = step.inputs.get('a') ?? 'undefined'; const b = step.inputs.get('b') ?? 'undefined';
    return `let ${step.variableName} = (${a}) > (${b});`;
  }
  case 'Greater Or Equal (>=)': {
    const a = step.inputs.get('a') ?? 'undefined'; const b = step.inputs.get('b') ?? 'undefined';
    return `let ${step.variableName} = (${a}) >= (${b});`;
  }
  case 'Less Than (<)': {
    const a = step.inputs.get('a') ?? 'undefined'; const b = step.inputs.get('b') ?? 'undefined';
    return `let ${step.variableName} = (${a}) < (${b});`;
  }
  case 'Less Or Equal (<=)': {
    const a = step.inputs.get('a') ?? 'undefined'; const b = step.inputs.get('b') ?? 'undefined';
    return `let ${step.variableName} = (${a}) <= (${b});`;
  }
  case 'Contains': {
    const t = step.inputs.get('text') ?? '\'\''; const s = step.inputs.get('substring') ?? '\'\'';
    return `let ${step.variableName} = (${t}).includes(${s});`;
  }
  case 'Starts With': {
    const t = step.inputs.get('text') ?? '\'\''; const p = step.inputs.get('prefix') ?? '\'\'';
    return `let ${step.variableName} = (${t}).startsWith(${p});`;
  }
  case 'Ends With': {
    const t = step.inputs.get('text') ?? '\'\''; const s = step.inputs.get('suffix') ?? '\'\'';
    return `let ${step.variableName} = (${t}).endsWith(${s});`;
  }
  case 'Is Null': {
    const v = step.inputs.get('value') ?? 'undefined';
    return `let ${step.variableName} = (${v}) == null;`;
  }
  case 'FromJSON': {
    const j = step.inputs.get('json') ?? '\'{}\'';
    return `let ${step.variableName} = JSON.parse(${j});`;
  }
  case 'ToJSON': {
    const v = step.inputs.get('value') ?? 'undefined';
    return `let ${step.variableName} = JSON.stringify(${v});`;
  }
  default:
    return `// Unknown utility: ${step.funcName}`;
  }
}

// ---------- instrumentation helpers ----------

function emitPreamble(runId: string): string[] {
  return [
    `const __ff_runId = '${runId}';`,
    `const __ff_ch = 'funcflow.exec.' + __ff_runId;`,
    'function __ff_summarize(v, declaredType) {',
    '  if (v == null) return {type:\'null\', value:null};',
    '  if (v.rowCount !== undefined && v.columns !== undefined)',
    '    return {type:\'dataframe\', rows:v.rowCount, cols:v.columns.length,',
    '      colNames:v.columns.names(), clone:v.clone()};',
    '  if (v.length !== undefined && v.name !== undefined && v.toList)',
    '    return {type:\'column\', name:v.name, length:v.length, sample:v.toList().slice(0,5)};',
    '  if (declaredType === \'graphics\' && typeof v === \'string\')',
    '    return {type:\'graphics\', value:v};',
    // A Viewer or Widget has a live DOM `.root`. The instrumented run is in the
    // same browser tab, so the event carries the object by reference (like the
    // DataFrame clone above) — keep it so the panel can mount `.root` directly.
    '  if (v.root != null && typeof Element !== \'undefined\' && v.root instanceof Element)',
    '    return {type: declaredType === \'viewer\' ? \'viewer\' : \'widget\', value:v};',
    '  if (typeof v === \'object\') return {type:\'object\', str:String(v).slice(0,200)};',
    '  return {type:\'primitive\', value:v};',
    '}',
    'function __ff_emit(type, nodeId, data) {',
    '  grok.events.fireCustomEvent(__ff_ch, Object.assign({type, nodeId, timestamp:Date.now()}, data||{}));',
    '}',
    `__ff_emit('run-start', '');`,
    '',
  ];
}

function wrapInstrumented(
  codeLine: string, step: CompiledStep, options: EmitOptions,
  extra?: {outputExpr?: string; declaredType?: string},
): string[] {
  const lines: string[] = [];
  lines.push(`__ff_emit('node-start', '${step.nodeId}');`);

  let bodyLine = codeLine;
  const letMatch = codeLine.match(/^let\s+(\w+)\s*=/);
  if (letMatch) {
    lines.push(`let ${letMatch[1]};`);
    bodyLine = codeLine.replace(/^let\s+/, '');
  }

  lines.push('try {');
  lines.push(`  ${bodyLine}`);
  if (extra?.outputExpr) {
    const typeArg = extra.declaredType ? `, '${extra.declaredType}'` : '';
    lines.push(`  __ff_emit('node-complete', '${step.nodeId}', ` +
      `{outputs:{${extra.outputExpr}: __ff_summarize(${extra.outputExpr}${typeArg})}});`);
  } else lines.push(`  __ff_emit('node-complete', '${step.nodeId}');`);
  lines.push('} catch (__ff_err) {');
  lines.push(`  __ff_emit('node-error', '${step.nodeId}', {error: __ff_err.message, stack: __ff_err.stack});`);
  if (options.haltOnError !== false) {
    lines.push(`  __ff_emit('run-complete', '', {success: false});`);
    lines.push('  throw __ff_err;');
  }
  lines.push('}');
  return lines;
}

function emitFuncStepInstrumented(step: CompiledStep, options: EmitOptions, flow: FlowEditor): string[] {
  const params: string[] = [];
  for (const [name, expr] of step.inputs) params.push(`${name}: ${expr}`);
  const paramsStr = params.length > 0 ? `{${params.join(', ')}}` : '{}';

  const lines: string[] = [];
  lines.push(`__ff_emit('node-start', '${step.nodeId}');`);
  lines.push(`let ${step.variableName};`);
  lines.push('try {');
  lines.push(`  ${step.variableName} = await grok.functions.call('${step.funcName}', ${paramsStr});`);

  // Build outputs summary.
  const node = flow.getNodeById(step.nodeId);
  const outputEntries: string[] = [];
  for (const [key, varName] of step.outputs) {
    const slotType = (node?.outputs as Record<string, {socket: {dgType: string}} | undefined>)[key]?.socket.dgType;
    const typeArg = slotType ? `, '${slotType}'` : '';
    outputEntries.push(`${varName}: __ff_summarize(${varName}${typeArg})`);
  }

  // In-place mutating function detection: dataframe input(s) but no real outputs.
  if (node && step.outputs.size === 0) {
    for (const [inKey, input] of Object.entries(node.inputs) as Array<[string, {socket: {dgType: string}} | undefined]>) {
      if (input?.socket.dgType !== 'dataframe') continue;
      const inputExpr = step.inputs.get(inKey);
      if (inputExpr && inputExpr !== 'undefined') {
        outputEntries.push(`'${inKey} (modified)': __ff_summarize(${inputExpr}, 'dataframe')`);
        break;
      }
    }
  }

  // SetVar produces no output, but its `value` input is the stored value —
  // surface it (labeled by the variable name) so clicking the node previews it
  // (table → grid, column → sample, …) exactly like any output-bearing node.
  if (node && node.dgFunc?.name?.toLowerCase() === 'setvar') {
    const valueExpr = step.inputs.get('value');
    if (valueExpr && valueExpr !== 'undefined') {
      const varLabel = JSON.stringify(String(node.inputValues['variableName'] ?? 'value'));
      const valueInput = (node.inputs as Record<string, {socket: {dgType: string}} | undefined>)['value'];
      const slotType = valueInput?.socket.dgType;
      const typeArg = slotType && slotType !== 'dynamic' ? `, '${slotType}'` : '';
      outputEntries.push(`${varLabel}: __ff_summarize(${valueExpr}${typeArg})`);
    }
  }

  const outputsObj = outputEntries.length > 0 ? `{${outputEntries.join(', ')}}` : '{}';
  lines.push(`  __ff_emit('node-complete', '${step.nodeId}', {outputs: ${outputsObj}});`);
  lines.push('} catch (__ff_err) {');
  lines.push(`  __ff_emit('node-error', '${step.nodeId}', {error: __ff_err.message, stack: __ff_err.stack});`);
  if (options.haltOnError !== false) {
    lines.push(`  __ff_emit('run-complete', '', {success: false});`);
    lines.push('  throw __ff_err;');
  }
  lines.push('}');
  return lines;
}

/** The `value` input expression of a SetVar func step, or null when the step
 *  isn't a SetVar or has no value. The slot type may be `dynamic` even when the
 *  value is a dataframe at runtime, so we can't decide at emit time whether the
 *  extra dataframe-name registration applies — the emitted code checks it at
 *  runtime (see emitSetVarByDataframeName). */
function setVarValueExpr(step: CompiledStep, flow: FlowEditor): string | null {
  const node = flow.getNodeById(step.nodeId);
  if (!node || node.dgFunc?.name?.toLowerCase() !== 'setvar') return null;
  const valueExpr = step.inputs.get('value');
  if (!valueExpr || valueExpr === 'undefined') return null;
  return valueExpr;
}

/** Runtime-guarded second `SetVar` registration: when the value is a dataframe,
 *  also register it keyed by the dataframe's runtime `.name`, so creation-script
 *  `GetVar`s that reference the table by its actual name resolve (the platform
 *  resolver registers both). The `instanceof` check is at runtime because the
 *  value slot is often `dynamic` even when it carries a dataframe. */
function emitSetVarByDataframeName(valueExpr: string): string {
  return `if (${valueExpr} instanceof DG.DataFrame) ` +
    `await grok.functions.call('SetVar', {variableName: ${valueExpr}.name, value: ${valueExpr}});`;
}

/** For each non-pass-through dataframe output, emit
 *  `if (varName != null) await varName.meta.detectSemanticTypes();` */
function emitDetectSemanticTypes(step: CompiledStep, flow: FlowEditor): string[] {
  const node = flow.getNodeById(step.nodeId);
  if (!node) return [];
  const lines: string[] = [];
  for (const [key, varName] of step.outputs) {
    if (key.endsWith(PASSTHROUGH_SUFFIX)) continue;
    const slotType = (node.outputs as Record<string, {socket: {dgType: string}} | undefined>)[key]?.socket.dgType;
    if (slotType === 'dataframe')
      lines.push(`if (${varName} != null) await ${varName}.meta.detectSemanticTypes();`);
  }
  return lines;
}

function emitBreakpointCode(step: CompiledStep): string[] {
  return [
    `__ff_emit('breakpoint-hit', '${step.nodeId}');`,
    'await new Promise((resolve) => {',
    `  const sub = grok.events.onCustomEvent(__ff_ch + '.continue').subscribe(() => {`,
    '    sub.unsubscribe();',
    '    resolve();',
    '  });',
    '});',
    `__ff_emit('node-complete', '${step.nodeId}');`,
  ];
}

function formatHeaderDefault(value: unknown, type: string): string {
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
