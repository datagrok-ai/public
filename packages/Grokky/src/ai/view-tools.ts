import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

export interface ViewToolDef {
  name: string;
  description: string;
  inputSchema?: object;
}

/** Tools sent to the Claude runtime: serializable defs go in the user message,
 * runners execute the matching input_request round-trips in the browser. */
export interface CollectedViewTools {
  defs: ViewToolDef[];
  runners: Map<string, (args: any) => any>;
}

export const NO_VIEW_TOOLS: CollectedViewTools = {defs: [], runners: new Map()};

const FIND_LIMIT = 10;

/** Functions applicable to a view, via `view.getFunctions()`: registered functions whose
 * `meta.viewType` matches, view-specific functions added by Dart overrides (query/script
 * editors), TableView commands, and — through the JsViewHost forwarding — functions a
 * JS-defined view (e.g. Flow) returns from its own `getFunctions()` override. */
function viewFunctions(view: DG.ViewBase | null): DG.Func[] {
  if (!view)
    return [];
  const funcs: DG.Func[] = [];
  const seen = new Set<string>();
  const collect = (list: any) => {
    if (Array.isArray(list))
      for (const f of list) {
        const name = f?.name;
        if (name && !seen.has(name)) {
          seen.add(name);
          funcs.push(f);
        }
      }
  };
  try {
    collect((view as any).getFunctions?.());
  } catch (e: any) {
    console.warn('Grokky: view.getFunctions failed:', e.message);
  }
  // Belt-and-braces for JS-defined views: reach the original ViewBase directly too.
  try {
    collect((view as any).jsView?.getFunctions?.());
  } catch (e: any) {
    console.warn('Grokky: jsView.getFunctions failed:', e.message);
  }
  return funcs;
}

/** OR-ranking: how many query words appear in the function's name/description.
 * AND-matching proved too strict — a query like "flow node pipeline graph" must
 * still surface listFlowNodes even though "pipeline" matches nothing. */
function matchScore(query: string, ...haystacks: (string | null | undefined)[]): number {
  const terms = query.toLowerCase().split(/\s+/).filter((t) => t.length > 0);
  if (terms.length === 0)
    return 1;
  const text = haystacks.filter(Boolean).join(' ').toLowerCase();
  return terms.filter((t) => text.includes(t)).length;
}

/** Compact, LLM-friendly description of a function: name, description, and inputs
 * (the `view` input is injected automatically at call time, so it is not listed). */
function funcBrief(f: DG.Func): object {
  let inputs: object[] = [];
  try {
    inputs = f.inputs
      .filter((p) => p.name !== 'view' && (p.propertyType as string) !== 'view')
      .map((p) => ({
        name: p.name,
        type: p.propertyType,
        ...(p.description ? {description: p.description} : {}),
      }));
  } catch (_) {}
  return {
    name: f.name,
    description: f.description || f.friendlyName || f.name,
    ...(inputs.length ? {inputs} : {}),
  };
}

/** Result values can be DG objects — reduce them to something JSON-serializable. */
function serializeResult(value: any): any {
  if (value == null || typeof value === 'string' || typeof value === 'number' || typeof value === 'boolean')
    return value;
  if (value instanceof DG.DataFrame)
    return {type: 'dataframe', name: value.name, rowCount: value.rowCount, columns: value.columns.names()};
  if (value instanceof DG.Column)
    return {type: 'column', name: value.name, length: value.length};
  if (value instanceof DG.ViewBase)
    return {type: 'view', name: value.name};
  try {
    return JSON.parse(JSON.stringify(value));
  } catch (_) {
    return String(value);
  }
}

async function invokeViewFunction(args: any): Promise<any> {
  const view = grok.shell.v;
  const name = args?.name;
  const funcs = viewFunctions(view);
  const f = funcs.find((x) => x.name === name);
  if (!f) {
    return {
      success: false,
      error: `No function '${name}' on the current view — call list_view_functions to see what is available`,
    };
  }
  const params: {[key: string]: any} = {...(args?.parameters ?? {})};
  for (const inp of f.inputs) {
    if (!(inp.name in params) && (inp.name === 'view' || (inp.propertyType as string) === 'view'))
      params[inp.name] = view;
  }
  const result = await f.apply(params);
  return {success: true, ...(result != null ? {result: serializeResult(result)} : {})};
}

const invocationSchema = {
  type: 'object',
  properties: {
    name: {type: 'string', description: 'Function name, exactly as returned by list_view_functions.'},
    parameters: {type: 'object', description: 'Arguments keyed by input name (the view argument is injected automatically).'},
  },
  required: ['name'],
};

/** The three meta-tools are static — same defs every turn (prompt-cache friendly);
 * the runners resolve the live current view at call time. */
const VIEW_TOOL_DEFS: ViewToolDef[] = [
  {
    name: 'list_view_functions',
    description: 'Search the functions applicable to the current view — its commands and view-specific ' +
      'operations (e.g. the query editor\'s SQL tools, Flow\'s graph tools). A view can have hundreds, ' +
      `so ALWAYS pass a query; at most ${FIND_LIMIT} matches are returned. ` +
      'Call this first when acting on the current view or when asked what you can do "here".',
    inputSchema: {
      type: 'object',
      properties: {
        query: {type: 'string', description: 'One or two broad words matched against function names and descriptions (OR-ranked — more matching words rank higher). Omit to see the first 10.'},
      },
    },
  },
  {
    name: 'get_view_function_result',
    description: 'Invoke a READ-ONLY function of the current view (one that inspects state without changing ' +
      'anything, e.g. get_query_info, list_flow_nodes) and return its result. ' +
      'For functions that change state, use call_view_function instead.',
    inputSchema: invocationSchema,
  },
  {
    name: 'call_view_function',
    description: 'Invoke a state-changing function of the current view (e.g. set_query_and_run, add_flow_node). ' +
      'Function name and parameters come from list_view_functions.',
    inputSchema: invocationSchema,
  },
];

/** Tool set declared to the runtime on every full-mode prompt. */
export function viewFunctionTools(): CollectedViewTools {
  const runners = new Map<string, (args: any) => any>();
  runners.set('list_view_functions', (args: any) => {
    const funcs = viewFunctions(grok.shell.v);
    let matched = funcs;
    if (args?.query) {
      matched = funcs
        .map((f) => ({f, score: matchScore(args.query, f.name, f.friendlyName, f.description)}))
        .filter((x) => x.score > 0)
        .sort((a, b) => b.score - a.score)
        .map((x) => x.f);
    }
    if (matched.length === 0 && funcs.length > 0)
      return {total: 0, functions: [], note: `No matches, but this view has ${funcs.length} functions — retry with a broader word or no query`};
    return {
      total: matched.length,
      functions: matched.slice(0, FIND_LIMIT).map(funcBrief),
      ...(matched.length > FIND_LIMIT ? {note: `Showing ${FIND_LIMIT} of ${matched.length} — refine the query`} : {}),
    };
  });
  runners.set('get_view_function_result', invokeViewFunction);
  runners.set('call_view_function', invokeViewFunction);
  return {defs: VIEW_TOOL_DEFS, runners};
}
