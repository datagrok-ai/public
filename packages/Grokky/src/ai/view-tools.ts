import * as DG from 'datagrok-api/dg';

export interface ViewToolDef {
  name: string;
  description: string;
  inputSchema?: object;
}

/** Tools collected for one prompt: serializable defs go to the Claude runtime,
 * runners execute the matching input_request round-trips in the browser. */
export interface CollectedViewTools {
  defs: ViewToolDef[];
  runners: Map<string, (args: any) => any>;
}

export const NO_VIEW_TOOLS: CollectedViewTools = {defs: [], runners: new Map()};

/** Collects the AI tools of a view: its own `getAITools()` (Dart-native tools, or a JS view's
 * override forwarded through the JsViewHost interop), plus tools from registered
 * `viewAIToolsProvider` functions whose `viewType` matches the view. */
export async function collectViewAITools(view: DG.ViewBase | null): Promise<CollectedViewTools> {
  if (!view)
    return NO_VIEW_TOOLS;
  const tools: DG.AIViewTool[] = [];
  try {
    const own = (view as any).getAITools?.();
    if (Array.isArray(own))
      tools.push(...own);
  } catch (e: any) {
    console.warn('Grokky: view.getAITools failed:', e.message);
  }
  for (const f of DG.Func.find({meta: {role: 'viewAIToolsProvider'}})) {
    if (f.options['viewType'] !== view.type)
      continue;
    try {
      const res = await f.apply({view});
      if (Array.isArray(res))
        tools.push(...res);
    } catch (e: any) {
      console.warn(`Grokky: viewAIToolsProvider ${f.name} failed:`, e.message);
    }
  }
  const defs: ViewToolDef[] = [];
  const runners = new Map<string, (args: any) => any>();
  for (const t of tools) {
    if (!t?.name || typeof t.run !== 'function' || runners.has(t.name))
      continue;
    runners.set(t.name, t.run);
    defs.push({name: t.name, description: t.description ?? t.name, ...(t.inputSchema ? {inputSchema: t.inputSchema} : {})});
  }
  return {defs, runners};
}
