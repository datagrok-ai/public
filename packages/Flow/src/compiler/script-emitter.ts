/** Generates the runnable Datagrok JavaScript from `CompiledStep[]`.
 *
 * Two emission modes:
 *   - clean: standard Datagrok script (with `//input:` / `//output:` headers)
 *   - instrumented: each step wrapped in try/catch firing `funcflow.exec.<runId>`
 *     events for live execution visualization. */

import {FlowEditor} from '../rete/flow-editor';
import {FlowNode, isExecKey, isSetVarNode} from '../rete/scheme';
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
  /** With `onlyNodeIds`, resolve connections whose source is *outside* the set
   *  to `_ffLive(...)` registry reads instead of in-script variables — so a
   *  single node can be re-run using values captured from a prior run. */
  liveExternalInputs?: boolean;
}

const PASSTHROUGH_SUFFIX = '__pt';

export function emitScript(
  flow: FlowEditor, settings: ScriptSettings, options?: EmitOptions,
): string {
  const liveBoundary = options?.liveExternalInputs && options.onlyNodeIds ? options.onlyNodeIds : undefined;
  const inst = options?.instrumented === true;
  // Instrumented runs clone every dataframe crossing into a func step — see
  // CompileOptions.cloneDataframeInputs (in-place functions must not mutate an
  // upstream node's captured value).
  let steps = compileGraph(flow, liveBoundary, {cloneDataframeInputs: inst});
  if (options?.onlyNodeIds) steps = steps.filter((s) => options.onlyNodeIds!.has(s.nodeId));
  const lines: string[] = [];

  lines.push(...buildHeaderLines(steps, flow, settings, 'javascript'));
  lines.push('');

  if (inst) lines.push(...emitPreamble(options!.runId!));

  // -------- body --------
  // Stash a step's live output values into the registry (instrumented only), so
  // a later single-node re-run can read them via `_ffLive`.
  const stash = (step: CompiledStep): void => {
    if (!inst) return;
    const line = stashLine(step, flow.getNodeById(step.nodeId));
    if (line) lines.push(line);
  };

  for (const step of steps) {
    if (step.nodeType === 'input') {
      stash(step); // param values are in body scope (declared by //input headers)
      continue;
    }

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
      // An Output node doubles as a SetVar: register the value in the run
      // context under the param name (and, for a dataframe, its runtime name)
      // so name-based consumers (Select Table, downstream scripts) resolve it —
      // Output and SetVar compile to the same thing.
      lines.push(`await grok.functions.call('SetVar', ` +
        `{variableName: ${JSON.stringify(step.variableName)}, value: ${inputExpr}});`);
      lines.push(emitSetVarByDataframeName(inputExpr));
      continue;
    }

    if (step.nodeType === 'utility') {
      if (step.funcName === 'Breakpoint') {
        if (inst && options?.enableBreakpoints) lines.push(...emitBreakpointCode(step));
        continue;
      }
      // Viewer nodes: `await table.plot.fromType(type, {}); v.setOptions(look)`.
      if (step.properties['viewerType']) {
        lines.push(...emitViewerStep(step, inst, options));
        stash(step);
        continue;
      }
      const code = emitUtilityStep(step);
      if (!code) continue;
      if (inst) {
        // Side-effect-only utilities (Log / Info / Warning) declare no variable
        // — summarizing `step.variableName` there would reference an undeclared
        // identifier and fail the node with a ReferenceError at run time. Only
        // steps that actually declare their variable get an output summary.
        const declaresVar = code.startsWith(`let ${step.variableName} =`);
        lines.push(...wrapInstrumented(code, step, options!,
          declaresVar ? {outputExpr: step.variableName} : undefined));
      } else
        lines.push(code);
      stash(step);
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
    if (setVarValue) {
      lines.push(emitSetVarByDataframeName(setVarValue));
      // …and a SetVar doubles as a script output: assign the stored value to
      // the output variable declared in the header (see setVarAsOutput) —
      // SetVar and Value Output compile to the same thing.
      const out = setVarAsOutput(step, flow);
      if (out) lines.push(`${out.name} = ${setVarValue};`);
    }

    // Auto-detect semantic types on dataframe outputs.
    lines.push(...emitDetectSemanticTypes(step, flow));
    stash(step);
  }

  if (inst) lines.push(`__ff_emit('run-complete', '', {success: true});`);

  return lines.join('\n');
}

/** The annotation-header block for the current steps: name / description /
 *  language / tags plus `//input:` / `//output:` lines derived from the
 *  Input/Output nodes (dataframe inputs first). Shared by the JS emission
 *  (`language: javascript`) and the `.flow` entity body (`language: flow`). */
function buildHeaderLines(
  steps: CompiledStep[], flow: FlowEditor, settings: ScriptSettings, language: string,
): string[] {
  const lines: string[] = [];
  lines.push(`//name: ${settings.name}`);
  if (settings.description) lines.push(`//description: ${settings.description}`);
  lines.push(`//language: ${language}`);
  if (settings.tags.length > 0) lines.push(`//tags: ${settings.tags.join(', ')}`);

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

  const emittedOutputs = new Set<string>();
  for (const step of steps.filter((s) => s.nodeType === 'output')) {
    const node = flow.getNodeById(step.nodeId);
    if (!node) continue;
    lines.push(buildOutputLine(step, node));
    emittedOutputs.add(String(step.properties['paramName']));
  }
  // SetVar nodes double as outputs — same `//output:` contract as a Value
  // Output wired to the same source, type inferred from the connection.
  for (const step of steps) {
    if (step.nodeType !== 'func') continue;
    const out = setVarAsOutput(step, flow);
    if (!out || emittedOutputs.has(out.name)) continue;
    emittedOutputs.add(out.name);
    lines.push(`//output: ${out.type} ${out.name}`);
  }
  return lines;
}

/** Header-only emission for a live editor — used by the `.flow` entity body,
 *  where the payload after the header is the ffjson document, not JS. */
export function emitHeaderLines(flow: FlowEditor, settings: ScriptSettings, language: string): string[] {
  return buildHeaderLines(compileGraph(flow), flow, settings, language);
}

/** `__ff_stash('<nodeId>', {<outputKey>: <valueExpr>, ...})` for a step's live
 *  outputs — real outputs (→ their variable) and dataframe passthroughs (→ the
 *  threaded, post-execution input value). Keyed by output socket key so a re-run
 *  can look up exactly the output a downstream connection reads. Null when the
 *  step produces nothing worth stashing (e.g. an Output node). */
function stashLine(step: CompiledStep, node: FlowNode | undefined): string | null {
  const entries: string[] = [];
  for (const [key, varName] of step.outputs)
    entries.push(`${JSON.stringify(key)}: ${varName}`);
  if (node) {
    for (const key of Object.keys(node.outputs)) {
      if (isExecKey(key) || !key.endsWith(PASSTHROUGH_SUFFIX)) continue;
      const inName = key.slice(0, -PASSTHROUGH_SUFFIX.length);
      const inExpr = step.inputs.get(inName);
      if (inExpr && inExpr !== 'undefined') entries.push(`${JSON.stringify(key)}: ${inExpr}`);
    }
  }
  if (entries.length === 0) return null;
  return `__ff_stash(${JSON.stringify(step.nodeId)}, {${entries.join(', ')}});`;
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

/** The named-argument literal for the `grok.functions.call` — a wrapper's
 *  reshaped arguments (`callInputs`) when present, else the node's inputs. */
function funcCallParamsStr(step: CompiledStep): string {
  const params: string[] = [];
  for (const [name, expr] of step.callInputs ?? step.inputs) params.push(`${name}: ${expr}`);
  return params.length > 0 ? `{${params.join(', ')}}` : '{}';
}

function emitFuncStep(step: CompiledStep): string {
  return `let ${step.variableName} = await grok.functions.call('${step.funcName}', ${funcCallParamsStr(step)});`;
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
    // Fail fast with the table name — a null table would otherwise surface as
    // a cryptic downstream error far from the node that caused it.
    const message = JSON.stringify(`Select Table: no open table or variable named "${raw}"`);
    return `let ${step.variableName} = ${expr};\n` +
      `if (${step.variableName} == null) throw new Error(${message});`;
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

// ---------- viewer nodes ----------

/** Non-empty look options, ready to pass to `setOptions` (drops blank values
 *  the user typed-then-cleared, and any leftover `#type` tag). */
function cleanViewerLook(look: unknown): Record<string, unknown> {
  const out: Record<string, unknown> = {};
  if (look && typeof look === 'object') {
    for (const [k, v] of Object.entries(look as Record<string, unknown>)) {
      if (k === '#type') continue;
      if (v === '' || v === null || v === undefined) continue;
      out[k] = v;
    }
  }
  return out;
}

/** The `setOptions(...)` object literal for a viewer step: each connected
 *  `column` option contributes `"<lookKey>": (<colExpr>).name` (a connected
 *  column wins over the panel value); the rest come from the stored look. */
function buildViewerOptions(step: CompiledStep): string {
  const specs = (step.properties['viewerOptionSpecs'] as Array<{key: string; kind: string}> | undefined) ?? [];
  const entries: string[] = [];
  const connected = new Set<string>();
  for (const s of specs) {
    if (s.kind !== 'column') continue;
    const expr = step.inputs.get(s.key);
    if (expr && expr !== 'undefined') {
      entries.push(`${JSON.stringify(s.key)}: (${expr}).name`);
      connected.add(s.key);
    }
  }
  for (const [k, val] of Object.entries(cleanViewerLook(step.properties['viewerLook']))) {
    if (connected.has(k)) continue; // a wired column overrides the typed value
    entries.push(`${JSON.stringify(k)}: ${JSON.stringify(val)}`);
  }
  return entries.length ? `{${entries.join(', ')}}` : '';
}

/** A viewer node → `let v = await <table>.plot.fromType('<Type>', {});` plus a
 *  `v.setOptions(<options>)` when options (typed or wired) are set. */
function emitViewerStep(step: CompiledStep, inst: boolean, options?: EmitOptions): string[] {
  const tableExpr = step.inputs.get('table');
  if (!tableExpr || tableExpr === 'undefined') return []; // no table wired — nothing to plot
  const type = String(step.properties['viewerType']);
  const optsLiteral = buildViewerOptions(step);
  const v = step.variableName;
  const create = `await ${tableExpr}.plot.fromType(${JSON.stringify(type)}, {})`;
  const setOpts = optsLiteral ? `${v}.setOptions(${optsLiteral});` : '';

  if (!inst) {
    const lines = [`let ${v} = ${create};`];
    if (setOpts) lines.push(setOpts);
    return lines;
  }

  const lines = [`__ff_emit('node-start', '${step.nodeId}');`, `let ${v};`, 'try {'];
  lines.push(`  ${v} = ${create};`);
  if (setOpts) lines.push(`  ${setOpts}`);
  // Slot-keyed like every other summary (see wrapInstrumented).
  const outKey = Array.from(step.outputs.entries()).find(([, varName]) => varName === v)?.[0] ?? v;
  lines.push(`  __ff_emit('node-complete', '${step.nodeId}', ` +
    `{outputs:{${JSON.stringify(outKey)}: __ff_summarize(${v}, 'viewer')}});`);
  lines.push('} catch (__ff_err) {');
  lines.push(`  __ff_emit('node-error', '${step.nodeId}', {error: __ff_err.message, stack: __ff_err.stack});`);
  if (options?.haltOnError !== false) {
    lines.push(`  __ff_emit('run-complete', '', {success: false});`);
    lines.push('  throw __ff_err;');
  }
  lines.push('}');
  return lines;
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
    '  if (v.length !== undefined && v.name !== undefined && v.toList) {',
    // Capture a one-column DataFrame (from a clone, so a later in-place mutation
    // of the source table can't change it) — the preview renders this as a grid
    // instead of a tiny text sample.
    '    var __ff_cdf = null;',
    '    try { __ff_cdf = DG.DataFrame.fromColumns([v.clone()]); } catch (e) {}',
    '    return {type:\'column\', name:v.name, length:v.length, sample:v.toList().slice(0,5), clone:__ff_cdf};',
    '  }',
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
    // Column-output summary with in-place detection: if the produced column is
    // (by instance) one of `inputTable`'s columns — the in-place idiom, e.g. Add
    // New Column mutates its input and returns the added column — capture a clone
    // of the whole table plus the column name, so the preview can show the column
    // IN the table it belongs to (scrolled into view) rather than alone.
    // Otherwise the base column summary (one-column grid / sample) is returned.
    'function __ff_col_summary(col, inputTable, declaredType) {',
    '  var base = __ff_summarize(col, declaredType);',
    '  if (base != null && base.type === \'column\' && inputTable != null &&',
    '      inputTable.columns !== undefined && inputTable.rowCount !== undefined) {',
    '    var owned = false;',
    '    try {',
    '      var __ff_cols = inputTable.columns.toList();',
    '      for (var __ff_i = 0; __ff_i < __ff_cols.length; __ff_i++) {',
    // Identity by the underlying Dart handle: each `.toList()`/`.col()` call
    // wraps the same column in a fresh JS object, so `===` on the wrappers can be
    // false even for the same column — compare `.dart`.
    '        var __ff_c = __ff_cols[__ff_i];',
    '        if (__ff_c === col || (__ff_c.dart != null && __ff_c.dart === col.dart)) { owned = true; break; }',
    '      }',
    '    } catch (e) {}',
    '    if (owned) {',
    '      try { base.tableClone = inputTable.clone(); base.scrollToColumn = col.name; } catch (e) {}',
    '    }',
    '  }',
    '  return base;',
    '}',
    'function __ff_emit(type, nodeId, data) {',
    '  grok.events.fireCustomEvent(__ff_ch, Object.assign({type, nodeId, timestamp:Date.now()}, data||{}));',
    '}',
    // Dims-only table summary for pass-through wires: enough for the on-edge
    // "N × K" count label, without the full clone __ff_summarize would take.
    'function __ff_dims(v) {',
    '  return (v != null && v.rowCount !== undefined && v.columns !== undefined) ?',
    '    {type:\'dataframe\', rows:v.rowCount, cols:v.columns.length} : null;',
    '}',
    // Snapshot a dataframe crossing into a step (anything else passes through).
    // Same duck-type check as __ff_summarize; a func step works on its own copy
    // so in-place transformations never mutate the upstream captured value.
    'function __ff_clone(v) {',
    '  return (v != null && v.rowCount !== undefined && v.columns !== undefined &&',
    '    typeof v.clone === \'function\') ? v.clone() : v;',
    '}',
    // Live-value registry (on the tab global): each node stashes its outputs so a
    // later single-node re-run can read them without re-running upstream.
    '  function __ff_stash(nodeId, map) {',
    '    (globalThis.__ffFlowLive = globalThis.__ffFlowLive || {})[nodeId] = map;',
    '  }',
    '  function _ffLive(nodeId, key) {',
    '    var vals = globalThis.__ffFlowLive && globalThis.__ffFlowLive[nodeId];',
    '    if (!vals || !(key in vals))',
    '      throw new Error(\'Flow: no captured value for \' + nodeId + \':\' + key + \' — run the flow first.\');',
    '    return vals[key];',
    '  }',
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
    // Key the summary by the output SLOT key, not the variable name — that's
    // what `labelOutgoingConnections` (edge counts) and the port-preview
    // lookup use (`state.outputs[outputKey]`).
    const slotKey = Array.from(step.outputs.entries())
      .find(([, varName]) => varName === extra.outputExpr)?.[0] ?? extra.outputExpr;
    lines.push(`  __ff_emit('node-complete', '${step.nodeId}', ` +
      `{outputs:{${JSON.stringify(slotKey)}: __ff_summarize(${extra.outputExpr}${typeArg})}});`);
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

/** The expression for the step's single connected dataframe input, or null when
 *  it has zero or more than one. Under `cloneDataframeInputs` this is the
 *  snapshot variable the call actually mutates — so a column output that lands
 *  in it is detectable by instance (see `__ff_col_summary`). */
function singleDataframeInputExpr(step: CompiledStep, node: FlowNode | undefined): string | null {
  if (!node) return null;
  const exprs: string[] = [];
  for (const [key, input] of Object.entries(node.inputs) as Array<[string, {socket: {dgType: string}} | undefined]>) {
    if (input?.socket.dgType !== 'dataframe') continue;
    const expr = step.inputs.get(key);
    if (expr && expr !== 'undefined') exprs.push(expr);
  }
  return exprs.length === 1 ? exprs[0] : null;
}

function emitFuncStepInstrumented(step: CompiledStep, options: EmitOptions, flow: FlowEditor): string[] {
  const paramsStr = funcCallParamsStr(step);

  const lines: string[] = [];
  lines.push(`__ff_emit('node-start', '${step.nodeId}');`);
  lines.push(`let ${step.variableName};`);
  // Dataframe input snapshots: declared in body scope (downstream steps read
  // the pass-through through them), assigned inside try (a `_ffLive` source
  // can throw — that must surface as THIS node's error).
  const snapshots = Array.from(step.cloneInputs ?? []);
  for (const [snapVar] of snapshots) lines.push(`let ${snapVar};`);
  lines.push('try {');
  for (const [snapVar, srcExpr] of snapshots) lines.push(`  ${snapVar} = __ff_clone(${srcExpr});`);
  lines.push(`  ${step.variableName} = await grok.functions.call('${step.funcName}', ${paramsStr});`);

  // Build outputs summary.
  const node = flow.getNodeById(step.nodeId);
  // The step's single connected dataframe input expression (the snapshot var
  // under cloneDataframeInputs) — the in-place table a column output may belong
  // to. Only when exactly one dataframe input is connected is the association
  // unambiguous; with zero or several we fall back to the plain column summary.
  const singleDfInput = singleDataframeInputExpr(step, node);
  const outputEntries: string[] = [];
  // Key each entry by the output *slot key* (what `labelOutgoingConnections`
  // and the live-stash look up), not by the value expression — a multi-output
  // func's expression is `<var>.<name>`, which is not a valid object-literal key.
  for (const [key, outExpr] of step.outputs) {
    const slotType = (node?.outputs as Record<string, {socket: {dgType: string}} | undefined>)[key]?.socket.dgType;
    const typeArg = slotType ? `, '${slotType}'` : '';
    if (slotType === 'column' && singleDfInput)
      outputEntries.push(`${JSON.stringify(key)}: __ff_col_summary(${outExpr}, ${singleDfInput}${typeArg})`);
    else
      outputEntries.push(`${JSON.stringify(key)}: __ff_summarize(${outExpr}${typeArg})`);
  }

  // Dataframe pass-throughs: dims-only summaries (rows × cols, no clone —
  // cheap) keyed by the pass-through slot key, so wires leaving a `<in> →`
  // port get their "N × K" count label just like wires from real outputs.
  // The panels skip `__pt` keys (value-inspector.ts).
  if (node) {
    for (const key of Object.keys(node.outputs)) {
      if (isExecKey(key) || !key.endsWith(PASSTHROUGH_SUFFIX)) continue;
      const slotType = (node.outputs as Record<string, {socket: {dgType: string}} | undefined>)[key]?.socket.dgType;
      if (slotType !== 'dataframe') continue;
      const inExpr = step.inputs.get(key.slice(0, -PASSTHROUGH_SUFFIX.length));
      if (inExpr && inExpr !== 'undefined')
        outputEntries.push(`${JSON.stringify(key)}: __ff_dims(${inExpr})`);
    }
  }

  // Capture the (possibly in-place-modified) table(s) threaded through a
  // dataframe input, keyed `<input> (modified)`. Covers two cases the real-output
  // summaries miss: a pure in-place mutator (no outputs at all), AND a node whose
  // real output isn't a table but still threads one through its passthrough —
  // e.g. AddNewColumn returns a *column*, yet a viewer wired to its "table →"
  // passthrough needs that post-execution table (for inspect + the column
  // picker). Every connected dataframe input is captured (a two-table mutator
  // with no output previews both), so the preview shows all it transformed.
  // Skipped when a real dataframe output already carries it.
  const hasDataframeOutput = Array.from(step.outputs.keys()).some((key) =>
    (node?.outputs as Record<string, {socket: {dgType: string}} | undefined>)[key]?.socket.dgType === 'dataframe');
  if (node && !hasDataframeOutput) {
    for (const [inKey, input] of Object.entries(node.inputs) as Array<[string, {socket: {dgType: string}} | undefined]>) {
      if (input?.socket.dgType !== 'dataframe') continue;
      const inputExpr = step.inputs.get(inKey);
      if (inputExpr && inputExpr !== 'undefined')
        outputEntries.push(`'${inKey} (modified)': __ff_summarize(${inputExpr}, 'dataframe')`);
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
  if (!node || !isSetVarNode(node)) return null;
  const valueExpr = step.inputs.get('value');
  if (!valueExpr || valueExpr === 'undefined') return null;
  return valueExpr;
}

const IDENTIFIER_RE = /^[A-Za-z_$][A-Za-z0-9_$]*$/;

/** The `//output:` contract a SetVar node contributes — SetVar and Value Output
 *  compile to the same thing, so a SetVar's variable becomes a script output of
 *  that name. Null when it can't be one: not a SetVar, the name is dynamic
 *  (wired in) or not a valid JS identifier (the run-context registration via
 *  the SetVar call still happens), or no value is connected. The type is
 *  inferred from the connected source socket — mirroring Value Output's
 *  on-connect auto-typing. */
function setVarAsOutput(step: CompiledStep, flow: FlowEditor): {name: string; type: string} | null {
  const node = flow.getNodeById(step.nodeId);
  if (!node || !isSetVarNode(node)) return null;
  const name = String(step.inputValues['variableName'] ?? '').trim();
  if (!IDENTIFIER_RE.test(name)) return null;
  const src = flow.getInputSource(step.nodeId, 'value');
  if (!src) return null;
  const dgType = (src.node.outputs as Record<string, {socket: {dgType: string}} | undefined>)[src.outputKey]
    ?.socket.dgType;
  return {name, type: dgType && dgType !== 'dynamic' && dgType !== 'object' ? dgType : 'dynamic'};
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
