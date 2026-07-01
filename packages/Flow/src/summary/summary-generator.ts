/** Heuristic plain-language summaries for nodes and whole flows (U12). A node's
 *  caption comes from, in order: a built-in type summary, a curated function
 *  summary, or a humanized fallback from the function's friendly name. The flow
 *  summary groups the graph into its disjoint pipelines (connected components)
 *  and describes each left-to-right. No platform calls — pure and testable. */

import {FlowNode, FlowConnection, isExecKey} from '../rete/scheme';
import {topologicalSortNodes} from '../compiler/topological-sort';
import {CURATED_FUNC_SUMMARIES, BUILTIN_SUMMARIES, str, prop} from './summary-defs';

/** Turn an identifier (camelCase / Pascal / snake / kebab / path) into a phrase:
 *  'addChemPropertiesColumns' → 'Add chem properties columns'; 'FilterRows' →
 *  'Filter rows'. */
export function humanize(name: string): string {
  const words = str(name)
    .replace(/[_\-./]+/g, ' ')
    .replace(/([a-z0-9])([A-Z])/g, '$1 $2')
    .replace(/([A-Z]+)([A-Z][a-z])/g, '$1 $2')
    .replace(/\s+/g, ' ')
    .trim();
  if (!words) return '';
  const lower = words.toLowerCase();
  return lower.charAt(0).toUpperCase() + lower.slice(1);
}

const clean = (s: string): string => str(s).replace(/\s{2,}/g, ' ');

/** The bare DG function name a node runs (strips any `Package:` prefix). */
function bareFuncName(node: FlowNode): string {
  return str(node.dgFunc?.name ?? node.dgFuncName).replace(/^.*:/, '');
}

/** A short, human caption for a single node. Never empty. */
export function summarizeNode(node: FlowNode): string {
  const t = str(node.dgTypeName);

  if (BUILTIN_SUMMARIES[t]) return clean(BUILTIN_SUMMARIES[t](node));
  if (t.startsWith('Inputs/')) {
    const p = prop(node, 'paramName') || str(node.label);
    const ty = str(node.dgOutputType);
    return clean(`Input “${p}”${ty ? ` (${ty})` : ''}`);
  }
  if (t.startsWith('Outputs/')) {
    const p = prop(node, 'paramName') || str(node.label);
    return clean(`Outputs “${p}”`);
  }
  if (t.startsWith('Constants/')) return clean(`Constant: ${str(node.label)}`);
  if (t.startsWith('Comparisons/')) return clean(`Compares: ${str(node.label)}`);

  const fname = bareFuncName(node);
  const curated = CURATED_FUNC_SUMMARIES[fname.toLowerCase()];
  if (curated) return clean(curated(node));

  // Fallback: a friendly name with spaces is already human; otherwise humanize.
  const friendly = str(node.dgFunc?.friendlyName);
  if (friendly && /\s/.test(friendly)) return clean(friendly);
  return clean(humanize(friendly || fname) || 'Runs a function');
}

/** Where one of a step's inputs comes from — another step in the same pipeline. */
export interface StepInput {
  /** The input port name (e.g. 'table', 'table1') — '' when unnamed. */
  key: string;
  /** 1-based index of the producing step within the pipeline. */
  from: number;
}
export interface FlowStep {
  /** 1-based position in the pipeline (dependency order). */
  index: number;
  caption: string;
  /** Upstream steps feeding this one — "what goes into where". */
  inputs: StepInput[];
}
export interface Pipeline {
  steps: FlowStep[];
}
export interface FlowSummary {
  /** One entry per disjoint pipeline (connected component). */
  pipelines: Pipeline[];
  /** A multi-line plain-text rendering (full, untruncated, dependency-ordered). */
  text: string;
  nodeCount: number;
  pipelineCount: number;
}

/** Summarize an entire flow. Steps are numbered in the **exact order the script
 *  executes** — the same `topologicalSortNodes` the compiler/emitter uses — and
 *  grouped into disjoint pipelines (which execution keeps contiguous, so each
 *  pipeline is a contiguous run of step numbers). Every step is annotated with
 *  which earlier step feeds each input ("what goes into where"). Captions are
 *  full — never truncated. */
export function summarizeFlow(nodes: FlowNode[], connections: FlowConnection[]): FlowSummary {
  if (nodes.length === 0)
    return {pipelines: [], text: 'Empty flow — nothing here yet.', nodeCount: 0, pipelineCount: 0};

  // The canonical execution order. On a cycle the sort returns a prefix, so
  // append any leftover nodes (original order) to keep every node numbered.
  const order = topologicalSortNodes(nodes, connections);
  const placed = new Set(order);
  const execOrder = [...order, ...nodes.filter((n) => !placed.has(n.id)).map((n) => n.id)];
  const stepNo = new Map(execOrder.map((id, i) => [id, i + 1])); // 1-based script line
  const byId = new Map(nodes.map((n) => [n.id, n]));

  // Disjoint pipelines (weakly connected components) via union-find.
  const parent = new Map<string, string>();
  for (const n of nodes) parent.set(n.id, n.id);
  const find = (x: string): string => {
    let r = x;
    while (parent.get(r) !== r) r = parent.get(r)!;
    while (parent.get(x) !== r) {
      const nx = parent.get(x)!;
      parent.set(x, r);
      x = nx;
    }
    return r;
  };
  for (const c of connections) {
    if (parent.has(c.source) && parent.has(c.target)) {
      const ra = find(c.source); const rb = find(c.target);
      if (ra !== rb) parent.set(ra, rb);
    }
  }
  const groups = new Map<string, string[]>();
  for (const id of execOrder) {
    const r = find(id);
    let arr = groups.get(r);
    if (!arr) {
      arr = [];
      groups.set(r, arr);
    }
    arr.push(id);
  }

  const inputsOf = (id: string): StepInput[] => connections
    .filter((c) => c.target === id && !isExecKey(str(c.targetInput)) && !isExecKey(str(c.sourceOutput)))
    .map((c) => ({key: str(c.targetInput), from: stepNo.get(c.source) ?? 0}))
    .filter((inp) => inp.from > 0 && inp.from !== stepNo.get(id))
    .sort((a, b) => a.from - b.from);

  // Pipelines ordered by their first executed step; steps within keep exec order.
  const pipelines: Pipeline[] = [...groups.values()]
    .map((ids) => ids.slice().sort((a, b) => stepNo.get(a)! - stepNo.get(b)!))
    .sort((a, b) => stepNo.get(a[0])! - stepNo.get(b[0])!)
    .map((ids) => ({
      steps: ids.map((id) => ({index: stepNo.get(id)!, caption: summarizeNode(byId.get(id)!), inputs: inputsOf(id)})),
    }));

  const multi = pipelines.length > 1;
  const renderInputs = (s: FlowStep): string => s.inputs.length === 0 ? '' :
    ` (← ${s.inputs.map((inp) => inp.key ? `${inp.key} from #${inp.from}` : `#${inp.from}`).join(', ')})`;
  const text = pipelines.map((p, k) => {
    const head = multi ? `Pipeline ${k + 1}:\n` : '';
    return head + p.steps.map((s) => `${multi ? '  ' : ''}${s.index}. ${s.caption}${renderInputs(s)}`).join('\n');
  }).join('\n\n');

  return {pipelines, text, nodeCount: nodes.length, pipelineCount: pipelines.length};
}
