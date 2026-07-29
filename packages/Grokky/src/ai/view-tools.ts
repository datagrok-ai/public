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
const WIDGET_LIMIT = 100;
const WIDGET_DEPTH = 4;

/** Defensive accessors — widgets can be Dart-backed wrappers or pure JS, old or new api. */
function widgetChildren(w: any): any[] {
  try {
    const c = w?.children;
    return Array.isArray(c) ? c.filter(Boolean) : [];
  } catch (_) {
    return [];
  }
}

function widgetAiDescription(w: any): string | null {
  try {
    const d = w?.aiDescription;
    if (d)
      return d;
  } catch (_) {}
  return null;
}

function widgetFunctions(w: any): DG.Func[] {
  try {
    const list = w?.getFunctions?.();
    return Array.isArray(list) ? list.filter((f: any) => f?.name) : [];
  } catch (_) {
    return [];
  }
}

export interface WidgetBrief {
  ref: string;
  type: string;
  aiDescription?: string;
  functions?: number;
}

/** Dialogs are widgets too, but they live outside the view's DOM subtree —
 * surface them as extra widget roots (refs 'dlg0', 'dlg0.1', ...). */
function openDialogs(): any[] {
  try {
    const list = (DG.Dialog as any).getOpenDialogs?.();
    return Array.isArray(list) ? list.filter(Boolean) : [];
  } catch (_) {
    return [];
  }
}

/** Flattens the widget tree of [view] plus any open dialogs into refs ('0', '0.2',
 * 'dlg0', 'dlg0.1', ...) — index paths through `Widget.children`, resolvable back
 * with resolveWidgetRef. */
export function collectWidgets(view: DG.ViewBase | null): WidgetBrief[] {
  const result: WidgetBrief[] = [];
  const walk = (w: any, ref: string, depth: number) => {
    if (depth > WIDGET_DEPTH || result.length >= WIDGET_LIMIT)
      return;
    if (ref !== '') {
      const ai = widgetAiDescription(w);
      const funcCount = widgetFunctions(w).length;
      result.push({
        ref,
        type: (() => { try { return w?.type ?? 'Widget'; } catch (_) { return 'Widget'; } })(),
        ...(ai ? {aiDescription: ai} : {}),
        ...(funcCount ? {functions: funcCount} : {}),
      });
    }
    widgetChildren(w).forEach((c, i) => walk(c, ref === '' ? `${i}` : `${ref}.${i}`, depth + 1));
  };
  if (view)
    walk(view, '', 0);
  openDialogs().forEach((d, i) => walk(d, `dlg${i}`, 1));
  return result;
}

/** Resolves a ref from collectWidgets back to the live widget; '' / null → the view itself. */
function resolveWidgetRef(view: DG.ViewBase | null, ref: string | null | undefined): any {
  if (ref == null || ref === '')
    return view;
  const parts = String(ref).split('.');
  const dlg = /^dlg(\d+)$/.exec(parts[0]);
  let w: any = view;
  if (dlg) {
    w = openDialogs()[Number(dlg[1])];
    parts.shift();
  }
  for (const part of parts) {
    if (!w)
      return null;
    w = widgetChildren(w)[Number(part)];
  }
  return w ?? null;
}

/** Widgets of [view] that carry an AI briefing — used in the workspace context snapshot. */
export function widgetBriefings(view: DG.ViewBase | null, max: number = 6): WidgetBrief[] {
  return collectWidgets(view).filter((w) => w.aiDescription).slice(0, max);
}

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

/** Function list for a tool call: the view's own set, or a sub-widget's when a
 * `widget` ref (from list_view_widgets) is passed. */
function targetFunctions(view: DG.ViewBase | null, widgetRef: string | null | undefined): DG.Func[] | null {
  if (widgetRef == null || widgetRef === '')
    return viewFunctions(view);
  const w = resolveWidgetRef(view, widgetRef);
  return w ? widgetFunctions(w) : null;
}

async function invokeViewFunction(args: any): Promise<any> {
  const view = grok.shell.v;
  const name = args?.name;
  const funcs = targetFunctions(view, args?.widget);
  if (funcs === null)
    return {success: false, error: `No widget '${args?.widget}' in the current view — call list_view_widgets for the current refs`};
  const f = funcs.find((x) => x.name === name);
  if (!f) {
    return {
      success: false,
      error: `No function '${name}' on the current ${args?.widget ? 'widget' : 'view'} — call list_view_functions to see what is available`,
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

const widgetRefProperty = {
  widget: {type: 'string', description: 'Optional widget ref from list_view_widgets (e.g. "0.2", or "dlg0" for an open dialog) to target that widget instead of the view itself.'},
};

const invocationSchema = {
  type: 'object',
  properties: {
    name: {type: 'string', description: 'Function name, exactly as returned by list_view_functions.'},
    parameters: {type: 'object', description: 'Arguments keyed by input name (the view argument is injected automatically).'},
    ...widgetRefProperty,
  },
  required: ['name'],
};

/** The meta-tools are static — same defs every turn (prompt-cache friendly);
 * the runners resolve the live current view (and optional widget ref) at call time. */
const VIEW_TOOL_DEFS: ViewToolDef[] = [
  {
    name: 'list_view_functions',
    description: 'Search the functions applicable to the current view — its commands and view-specific ' +
      'operations (e.g. the query editor\'s SQL tools, Flow\'s graph tools). A view can have hundreds, ' +
      `so ALWAYS pass a query; at most ${FIND_LIMIT} matches are returned. ` +
      'Call this first when acting on the current view or when asked what you can do "here". ' +
      'Pass a widget ref (from list_view_widgets) to list a specific sub-widget\'s functions instead.',
    inputSchema: {
      type: 'object',
      properties: {
        query: {type: 'string', description: 'One or two broad words matched against function names and descriptions (OR-ranked — more matching words rank higher). Omit to see the first 10.'},
        ...widgetRefProperty,
      },
    },
  },
  {
    name: 'list_view_widgets',
    description: 'List the widget tree of the current view (grid, viewers, editors, tab controls, panels) ' +
      'plus any OPEN DIALOGS (refs dlg0, dlg0.1, ...), each with its ref, type, AI briefing ' +
      '(aiDescription), and function count. Use a ref as the `widget` argument of list_view_functions / ' +
      'get_view_function_result / call_view_function to inspect or drive that specific widget — e.g. a ' +
      'dialog exposes get_dialog_info / set_input / click_button. Call it when the view-level functions ' +
      'are not enough, or whenever a dialog is open.',
    inputSchema: {type: 'object', properties: {}},
  },
  {
    name: 'get_view_function_result',
    description: 'Invoke a READ-ONLY function of the current view or one of its widgets (one that inspects ' +
      'state without changing anything, e.g. get_query_info, list_flow_nodes) and return its result. ' +
      'For functions that change state, use call_view_function instead.',
    inputSchema: invocationSchema,
  },
  {
    name: 'call_view_function',
    description: 'Invoke a state-changing function of the current view or one of its widgets ' +
      '(e.g. set_query_and_run, add_flow_node). Function name and parameters come from list_view_functions.',
    inputSchema: invocationSchema,
  },
];

/** Tool set declared to the runtime on every full-mode prompt. */
export function viewFunctionTools(): CollectedViewTools {
  const runners = new Map<string, (args: any) => any>();
  runners.set('list_view_widgets', () => {
    const widgets = collectWidgets(grok.shell.v);
    return {
      total: widgets.length,
      widgets,
      ...(widgets.length === 0 ? {note: 'This view exposes no sub-widgets — use the view-level functions'} : {}),
    };
  });
  runners.set('list_view_functions', (args: any) => {
    const funcs = targetFunctions(grok.shell.v, args?.widget);
    if (funcs === null)
      return {total: 0, functions: [], note: `No widget '${args?.widget}' in the current view — call list_view_widgets for the current refs`};
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
