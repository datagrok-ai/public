import type {FlowEditor} from './rete/flow-editor';
import type {FlowNode} from './rete/scheme';
import {missingRequiredInputs} from './rete/scheme';
import type {ExecutionController} from './execution/execution-controller';
import {NodeExecStatus} from './execution/execution-state';
import {ensureFunctionsRegistered, getRegisteredFuncs, getRegisteredTypeNames,
  getInputTypesForType, getOutputTypesForType,
  findNodeTypesAcceptingInput, findNodeTypesProducingOutput} from './rete/node-factory';
import {areTypesCompatible} from './types/type-map';
import {validateGraph} from './compiler/validator';
import {inputBlockReason} from './utils/input-values';
import type {Guide} from './guide/guide-model';
import {TUTORIALS, QUESTIONS} from './guide/guide-content';

/** Facade over FuncFlowView privates handed to the AI view functions (package.ts).
 *  The functions receive the generic current view and reach the FuncFlowView
 *  instance through `view.jsView.aiContext()`. */
export interface FlowAIContext {
  flow(): FlowEditor;
  execution(): ExecutionController | null;
  addNodeByType(typeName: string): Promise<FlowNode | null>;
  run(): void;
  /** Fire-and-forget: guides wait for user actions, so the function must not await them. */
  runGuide(guide: Guide): void;
}

const RUN_TIMEOUT_MS = 100_000;
const DEFAULT_FIND_LIMIT = 15;

/** Resolves the FuncFlowView behind the generic view wrapper (grok.shell.v). */
function flowCtx(view: any): FlowAIContext {
  const v = view?.jsView ?? view;
  const ctx: FlowAIContext | null = typeof v?.aiContext === 'function' ? v.aiContext() : null;
  if (!ctx)
    throw new Error('The current view is not a Flow editor');
  return ctx;
}

function nodeRef(ctx: FlowAIContext, id: string): FlowNode {
  const node = ctx.flow().getNodeById(id);
  if (!node)
    throw new Error(`Node '${id}' not found — call listFlowNodes for current ids`);
  return node;
}

function portMap(slots: Record<string, any>): Record<string, string> {
  const res: Record<string, string> = {};
  for (const [key, slot] of Object.entries<any>(slots ?? {}))
    res[key] = slot?.socket?.dgType ?? 'dynamic';
  return res;
}

function nodeBrief(flow: FlowEditor, exec: ExecutionController | null, n: FlowNode) {
  const state = exec?.state.getNodeState(n.id);
  const values = Object.entries(n.inputValues ?? {}).filter(([, v]) => v !== '' && v != null);
  return {
    id: n.id,
    label: n.label,
    typeName: n.dgTypeName,
    ...(n.dgFuncName ? {func: n.dgFuncName} : {}),
    ...(state ? {status: state.status, ...(state.error ? {error: state.error} : {})} : {}),
    ...(values.length ? {inputValues: Object.fromEntries(values)} : {}),
  };
}

function summarizeOutputs(outputs?: Record<string, any>): Record<string, any> | undefined {
  if (!outputs)
    return undefined;
  const res: Record<string, any> = {};
  for (const [key, s] of Object.entries<any>(outputs)) {
    res[key] = s?.type === 'dataframe' ? {type: 'dataframe', rows: s.rows, cols: s.cols} :
      s?.type === 'column' ? {type: 'column', name: s.name, length: s.length} :
      s?.type === 'scalar' ? {type: 'scalar', value: s.value} :
      {type: s?.type ?? 'unknown'};
  }
  return res;
}

function matchesQuery(query: string, ...haystacks: (string | undefined)[]): boolean {
  const terms = query.toLowerCase().split(/\s+/).filter((t) => t.length > 0);
  const text = haystacks.filter(Boolean).join(' ').toLowerCase();
  return terms.every((t) => text.includes(t));
}

// --- implementations of the registered Flow view functions (see package.ts) ---

export function listFlowNodes(view: any) {
  const ctx = flowCtx(view);
  const flow = ctx.flow();
  const exec = ctx.execution();
  return {
    nodes: flow.getNodes().map((n) => nodeBrief(flow, exec, n)),
    connections: flow.getConnections().map((c) => `${c.source}.${c.sourceOutput} -> ${c.target}.${c.targetInput}`),
  };
}

export function getFlowNodeDetails(view: any, nodeId: string) {
  const ctx = flowCtx(view);
  const flow = ctx.flow();
  const node = nodeRef(ctx, nodeId);
  const state = ctx.execution()?.state.getNodeState(node.id);
  const connected = (key: string) => flow.isInputConnected(node.id, key);
  return {
    id: node.id,
    label: node.label,
    typeName: node.dgTypeName,
    ...(node.dgFuncName ? {func: node.dgFuncName} : {}),
    inputs: portMap((node as any).inputs),
    outputs: portMap((node as any).outputs), // keys ending in __pt are pass-throughs of the same-named input
    inputValues: node.inputValues ?? {},
    ...(node.dgNodeType === 'input' ? {
      // The configured parameter value (set on the node / context panel);
      // when unresolvable, `valueBlocks` says why runs would need a dialog.
      value: node.properties['defaultValue'] ?? '',
      ...(inputBlockReason(node) ? {valueBlocks: inputBlockReason(node)} : {}),
    } : {}),
    missingRequired: missingRequiredInputs(node, connected),
    ...(state ? {
      status: state.status,
      ...(state.error ? {error: state.error} : {}),
      ...(summarizeOutputs(state.outputs) ? {lastOutputs: summarizeOutputs(state.outputs)} : {}),
    } : {}),
  };
}

export function findFlowNodeTypes(view: any, query?: string, acceptsInputType?: string,
  producesOutputType?: string, limit?: number) {
  flowCtx(view);
  ensureFunctionsRegistered();
  const max = limit ?? DEFAULT_FIND_LIMIT;
  let candidates: {typeName: string; label: string; description?: string}[];
  if (acceptsInputType)
    candidates = findNodeTypesAcceptingInput(acceptsInputType);
  else if (producesOutputType)
    candidates = findNodeTypesProducingOutput(producesOutputType);
  else {
    const funcTypes = getRegisteredFuncs().map((f) => ({
      typeName: f.nodeTypeName,
      label: f.name,
      description: [f.category, f.func?.description].filter(Boolean).join(' — '),
    }));
    const funcTypeNames = new Set(funcTypes.map((f) => f.typeName));
    const builtins = getRegisteredTypeNames()
      .filter((t) => !funcTypeNames.has(t))
      .map((t) => ({typeName: t, label: t.split('/').pop() ?? t, description: undefined as string | undefined}));
    candidates = [...funcTypes, ...builtins];
  }
  if (query)
    candidates = candidates.filter((c) => matchesQuery(query, c.typeName, c.label, c.description));
  if (producesOutputType && acceptsInputType) {
    candidates = candidates.filter((c) => {
      const outs = getOutputTypesForType(c.typeName);
      return [...outs.real, ...outs.passthrough].some((t) => areTypesCompatible(t, producesOutputType));
    });
  }
  return {
    total: candidates.length,
    matches: candidates.slice(0, max).map((c) => ({
      typeName: c.typeName,
      label: c.label,
      ...(c.description ? {description: c.description} : {}),
      inputTypes: getInputTypesForType(c.typeName),
      outputTypes: getOutputTypesForType(c.typeName).real,
    })),
  };
}

export async function addFlowNode(view: any, typeName: string, label?: string, inputValues?: object) {
  const ctx = flowCtx(view);
  const node = await ctx.addNodeByType(typeName);
  if (!node)
    return {success: false, error: `Unknown node type: ${typeName}`};
  if (label)
    node.label = label;
  if (inputValues && typeof inputValues === 'object') {
    Object.assign(node.inputValues, inputValues);
    ctx.flow().notifyNodeParamsChanged(node.id);
  }
  await ctx.flow().updateNode(node.id);
  return {
    success: true,
    nodeId: node.id,
    inputs: portMap((node as any).inputs),
    outputs: portMap((node as any).outputs),
  };
}

export async function connectFlowNodes(view: any, sourceNodeId: string, sourceOutput: string,
  targetNodeId: string, targetInput: string) {
  const ctx = flowCtx(view);
  const flow = ctx.flow();
  const source = nodeRef(ctx, sourceNodeId);
  const target = nodeRef(ctx, targetNodeId);
  const out = (source as any).outputs?.[sourceOutput];
  const inp = (target as any).inputs?.[targetInput];
  if (!out)
    return {success: false, error: `Output '${sourceOutput}' not found; available: ${Object.keys((source as any).outputs ?? {}).join(', ')}`};
  if (!inp)
    return {success: false, error: `Input '${targetInput}' not found; available: ${Object.keys((target as any).inputs ?? {}).join(', ')}`};
  const srcType = out.socket?.dgType ?? 'dynamic';
  const tgtType = inp.socket?.dgType ?? 'dynamic';
  if (!areTypesCompatible(srcType, tgtType))
    return {success: false, error: `Incompatible types: ${srcType} -> ${tgtType}`};
  await flow.addConnectionByKeys(sourceNodeId, sourceOutput, targetNodeId, targetInput);
  return {success: true};
}

export async function setFlowNodeInputs(view: any, nodeId: string, values: object) {
  const ctx = flowCtx(view);
  const node = nodeRef(ctx, nodeId);
  // On an input node, `value` means the configured parameter value (the same
  // store the on-node editor writes) — with it set, runs need no dialog.
  const vals: Record<string, unknown> = {...(values ?? {})};
  if (node.dgNodeType === 'input' && 'value' in vals) {
    node.properties['defaultValue'] = vals['value'];
    delete vals['value'];
  }
  Object.assign(node.inputValues, vals);
  ctx.flow().notifyNodeParamsChanged(node.id);
  await ctx.flow().updateNode(node.id);
  return {success: true, inputValues: node.inputValues,
    ...(node.dgNodeType === 'input' ? {value: node.properties['defaultValue'] ?? ''} : {})};
}

export async function selectFlowNode(view: any, nodeId: string) {
  const ctx = flowCtx(view);
  nodeRef(ctx, nodeId);
  await ctx.flow().selectNode(nodeId);
  return {success: true};
}

export function listFlowGuides(view: any, query?: string) {
  flowCtx(view);
  let guides = [...TUTORIALS, ...QUESTIONS];
  if (query)
    guides = guides.filter((g) => matchesQuery(query, g.id, g.title, g.summary));
  return guides.map((g) => ({id: g.id, kind: g.kind, title: g.title, ...(g.summary ? {summary: g.summary} : {})}));
}

export function startFlowGuide(view: any, guideId: string) {
  const ctx = flowCtx(view);
  const guide = [...TUTORIALS, ...QUESTIONS].find((g) => g.id === guideId);
  if (!guide)
    return {success: false, error: `Unknown guide '${guideId}' — call listFlowGuides for ids`};
  ctx.runGuide(guide);
  return {success: true, started: guide.title};
}

export async function runFlow(view: any) {
  const ctx = flowCtx(view);
  const flow = ctx.flow();
  const exec = ctx.execution();
  if (!exec)
    return {success: false, error: 'Execution is not available in this view'};
  const problems = validateGraph(flow);
  if (problems.some((p) => p.severity === 'error'))
    return {success: false, validation: problems};
  const done = new Promise<boolean>((resolve) => {
    const prev = exec.onRunEnd;
    exec.onRunEnd = (success: boolean) => {
      exec.onRunEnd = prev;
      prev?.(success);
      resolve(success);
    };
  });
  ctx.run();
  const timeout = new Promise<null>((resolve) => setTimeout(() => resolve(null), RUN_TIMEOUT_MS));
  const success = await Promise.race([done, timeout]);
  if (success === null)
    return {success: false, error: `Run did not finish within ${RUN_TIMEOUT_MS / 1000}s (may still be running)`};
  const failed = flow.getNodes()
    .map((n) => ({n, state: exec.state.getNodeState(n.id)}))
    .filter(({state}) => state?.status === NodeExecStatus.errored)
    .map(({n, state}) => ({id: n.id, label: n.label, error: state!.error}));
  return {success, ...(failed.length ? {failedNodes: failed} : {})};
}
