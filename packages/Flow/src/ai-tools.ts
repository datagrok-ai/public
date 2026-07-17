import * as DG from 'datagrok-api/dg';
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
import type {Guide} from './guide/guide-model';
import {TUTORIALS, QUESTIONS} from './guide/guide-content';

/** Facade over FuncFlowView privates handed to the tool builder. */
export interface FlowAIContext {
  flow(): FlowEditor;
  execution(): ExecutionController | null;
  addNodeByType(typeName: string): Promise<FlowNode | null>;
  run(): void;
  /** Fire-and-forget: guides wait for user actions, so the tool must not await them. */
  runGuide(guide: Guide): void;
}

const RUN_TIMEOUT_MS = 100_000;
const DEFAULT_FIND_LIMIT = 15;

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

export function buildFlowAITools(ctx: FlowAIContext): DG.AIViewTool[] {
  const nodeRef = (id: string): FlowNode => {
    const node = ctx.flow().getNodeById(id);
    if (!node)
      throw new Error(`Node '${id}' not found — call list_flow_nodes for current ids`);
    return node;
  };

  return [
    {
      name: 'list_flow_nodes',
      description: 'List the current flow graph: all nodes (id, label, type, status, set input values) ' +
        'and connections. Call this first to understand what is on the canvas.',
      run: () => {
        const flow = ctx.flow();
        const exec = ctx.execution();
        return {
          nodes: flow.getNodes().map((n) => nodeBrief(flow, exec, n)),
          connections: flow.getConnections().map((c) => `${c.source}.${c.sourceOutput} -> ${c.target}.${c.targetInput}`),
        };
      },
    },
    {
      name: 'get_flow_node_details',
      description: 'Ports (with DG types), editable input values, unmet requirements, and last-run ' +
        'outputs of one node.',
      inputSchema: {type: 'object', properties: {nodeId: {type: 'string'}}, required: ['nodeId']},
      run: (args: any) => {
        const flow = ctx.flow();
        const node = nodeRef(args.nodeId);
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
          missingRequired: missingRequiredInputs(node, connected),
          ...(state ? {
            status: state.status,
            ...(state.error ? {error: state.error} : {}),
            ...(summarizeOutputs(state.outputs) ? {lastOutputs: summarizeOutputs(state.outputs)} : {}),
          } : {}),
        };
      },
    },
    {
      name: 'find_flow_node_types',
      description: 'Search the flow node catalog (a curated subset of platform functions plus ' +
        'input/output/utility nodes). ALWAYS filter: pass `query` with what the node should do ' +
        '(e.g. "join tables", "open file", "mutate molecules"), and/or a DG type it must accept ' +
        'or produce (dataframe, column, string, ...). Returns at most `limit` (default 15) matches ' +
        'with their input/output types.',
      inputSchema: {
        type: 'object',
        properties: {
          query: {type: 'string', description: 'Words describing what the node does.'},
          acceptsInputType: {type: 'string', description: 'DG type one of its inputs must accept.'},
          producesOutputType: {type: 'string', description: 'DG type one of its outputs must produce.'},
          limit: {type: 'integer'},
        },
      },
      run: (args: any) => {
        ensureFunctionsRegistered();
        const limit = args?.limit ?? DEFAULT_FIND_LIMIT;
        let candidates: {typeName: string; label: string; description?: string}[];
        if (args?.acceptsInputType)
          candidates = findNodeTypesAcceptingInput(args.acceptsInputType);
        else if (args?.producesOutputType)
          candidates = findNodeTypesProducingOutput(args.producesOutputType);
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
        if (args?.query)
          candidates = candidates.filter((c) => matchesQuery(args.query, c.typeName, c.label, c.description));
        if (args?.producesOutputType && args?.acceptsInputType) {
          candidates = candidates.filter((c) => {
            const outs = getOutputTypesForType(c.typeName);
            return [...outs.real, ...outs.passthrough].some((t) => areTypesCompatible(t, args.producesOutputType));
          });
        }
        return {
          total: candidates.length,
          matches: candidates.slice(0, limit).map((c) => ({
            typeName: c.typeName,
            label: c.label,
            ...(c.description ? {description: c.description} : {}),
            inputTypes: getInputTypesForType(c.typeName),
            outputTypes: getOutputTypesForType(c.typeName).real,
          })),
        };
      },
    },
    {
      name: 'add_flow_node',
      description: 'Add a node to the canvas by its registered typeName (from find_flow_node_types). ' +
        'Optionally set editable input values right away. Returns the new node id and its ports.',
      inputSchema: {
        type: 'object',
        properties: {
          typeName: {type: 'string'},
          label: {type: 'string', description: 'Optional custom title.'},
          inputValues: {type: 'object', description: 'Editable primitive inputs, key → value.'},
        },
        required: ['typeName'],
      },
      run: async (args: any) => {
        const node = await ctx.addNodeByType(args.typeName);
        if (!node)
          return {success: false, error: `Unknown node type: ${args.typeName}`};
        if (args.label)
          node.label = args.label;
        if (args.inputValues && typeof args.inputValues === 'object') {
          Object.assign(node.inputValues, args.inputValues);
          ctx.flow().notifyNodeParamsChanged(node.id);
        }
        await ctx.flow().updateNode(node.id);
        return {
          success: true,
          nodeId: node.id,
          inputs: portMap((node as any).inputs),
          outputs: portMap((node as any).outputs),
        };
      },
    },
    {
      name: 'connect_flow_nodes',
      description: 'Connect a source node output to a target node input (port keys from ' +
        'get_flow_node_details / add_flow_node). Types must be compatible.',
      inputSchema: {
        type: 'object',
        properties: {
          sourceNodeId: {type: 'string'},
          sourceOutput: {type: 'string'},
          targetNodeId: {type: 'string'},
          targetInput: {type: 'string'},
        },
        required: ['sourceNodeId', 'sourceOutput', 'targetNodeId', 'targetInput'],
      },
      run: async (args: any) => {
        const flow = ctx.flow();
        const source = nodeRef(args.sourceNodeId);
        const target = nodeRef(args.targetNodeId);
        const out = (source as any).outputs?.[args.sourceOutput];
        const inp = (target as any).inputs?.[args.targetInput];
        if (!out)
          return {success: false, error: `Output '${args.sourceOutput}' not found; available: ${Object.keys((source as any).outputs ?? {}).join(', ')}`};
        if (!inp)
          return {success: false, error: `Input '${args.targetInput}' not found; available: ${Object.keys((target as any).inputs ?? {}).join(', ')}`};
        const srcType = out.socket?.dgType ?? 'dynamic';
        const tgtType = inp.socket?.dgType ?? 'dynamic';
        if (!areTypesCompatible(srcType, tgtType))
          return {success: false, error: `Incompatible types: ${srcType} -> ${tgtType}`};
        await flow.addConnectionByKeys(args.sourceNodeId, args.sourceOutput, args.targetNodeId, args.targetInput);
        return {success: true};
      },
    },
    {
      name: 'set_flow_node_inputs',
      description: 'Set editable input values of a node (key → value; keys from get_flow_node_details ' +
        'inputValues/inputs). Marks the node and its downstream stale.',
      inputSchema: {
        type: 'object',
        properties: {nodeId: {type: 'string'}, values: {type: 'object'}},
        required: ['nodeId', 'values'],
      },
      run: async (args: any) => {
        const node = nodeRef(args.nodeId);
        Object.assign(node.inputValues, args.values ?? {});
        ctx.flow().notifyNodeParamsChanged(node.id);
        await ctx.flow().updateNode(node.id);
        return {success: true, inputValues: node.inputValues};
      },
    },
    {
      name: 'select_flow_node',
      description: 'Select a node on the canvas so the user sees it (opens its properties panel).',
      inputSchema: {type: 'object', properties: {nodeId: {type: 'string'}}, required: ['nodeId']},
      run: async (args: any) => {
        nodeRef(args.nodeId);
        await ctx.flow().selectNode(args.nodeId);
        return {success: true};
      },
    },
    {
      name: 'list_flow_guides',
      description: 'List Flow\'s built-in interactive guides: step-by-step tutorials and short ' +
        '"how do I…" walkthroughs that highlight the actual UI. When the user asks how to do ' +
        'something in Flow, check here first — a matching guide beats a textual explanation.',
      inputSchema: {type: 'object', properties: {query: {type: 'string', description: 'Words to filter by.'}}},
      run: (args: any) => {
        let guides = [...TUTORIALS, ...QUESTIONS];
        if (args?.query)
          guides = guides.filter((g) => matchesQuery(args.query, g.id, g.title, g.summary));
        return guides.map((g) => ({id: g.id, kind: g.kind, title: g.title, ...(g.summary ? {summary: g.summary} : {})}));
      },
    },
    {
      name: 'start_flow_guide',
      description: 'Start an interactive guide (id from list_flow_guides) — it highlights the real ' +
        'UI step by step and waits for the user to act. ALWAYS confirm with the user first ' +
        '(AskUserQuestion) before starting one; never launch it unasked.',
      inputSchema: {type: 'object', properties: {guideId: {type: 'string'}}, required: ['guideId']},
      run: (args: any) => {
        const guide = [...TUTORIALS, ...QUESTIONS].find((g) => g.id === args.guideId);
        if (!guide)
          return {success: false, error: `Unknown guide '${args.guideId}' — call list_flow_guides for ids`};
        ctx.runGuide(guide);
        return {success: true, started: guide.title};
      },
    },
    {
      name: 'run_flow',
      description: 'Validate and execute the whole flow. Returns validation problems instead of running ' +
        'if the graph is invalid; otherwise waits for the run and reports per-node failures.',
      run: async () => {
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
      },
    },
  ];
}
