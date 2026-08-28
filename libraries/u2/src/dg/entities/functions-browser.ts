/* The platform feed for the core FunctionsBrowser (the Dart FunctionsWidget as a u2 factory):
   entries off the synchronous client registry (`DG.Func.find({})`), the roles pane off
   `DG.functionRoles`, the platform icon through the object handler, run through
   `prepare().edit()`, and selection into `grok.shell.o`. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {FunctionsBrowser, FunctionsBrowserOptions, FuncItem}
  from '../../components/collections/functions-browser.js';

export interface DgFunctionsBrowserOptions extends
  Omit<FunctionsBrowserOptions, 'items' | 'roles' | 'runAction' | 'icon'> {
  filter?: (f: DG.Func) => boolean;
  /** Only functions with one scalar (or dynamic) output — the Dart widget mode. */
  scalarOnly?: boolean;
  /** Packages whose functions are left out (`namespace.startsWith('Pkg:')` in Dart). */
  ignorePackages?: string[];
  /** Selection makes the function the shell's current object; default true. */
  setCurrentObject?: boolean;
}

/** `Types.isScalar` (ddt utils.dart:180) plus `dynamic`, as the Dart widget-mode predicate. */
const SCALAR_TYPES = new Set(['int', 'float', 'double', 'num', 'qnum', 'bigint',
  'bool', 'string', 'datetime', 'dynamic']);

function hasScalarOutput(f: DG.Func): boolean {
  const outputs = f.outputs;
  return outputs.length === 1 && SCALAR_TYPES.has(outputs[0].propertyType);
}

/** `(in, puts) : outType` — the Dart `renderParams` string. */
function signatureOf(f: DG.Func): string {
  const inputs = f.inputs.map((p) => p.name).join(', ');
  const outputs = f.outputs.map((p) => p.propertyType).join(', ');
  return `(${inputs}) : ${outputs === '' ? 'void' : outputs}`;
}

function funcItem(f: DG.Func): FuncItem {
  const name = f.nqName;
  const label = f.friendlyName || f.name;
  const description = f.description ?? '';
  const tags = (f.tags ?? []).map((t) => t.toLowerCase());
  // the Dart role set (`hasTagOrRole`): the tags plus the comma-separated `role` option
  const role = f.options['role'] as string | null | undefined;
  const roles = !role ? tags :
    [...tags, ...role.split(',').map((r) => r.trim().toLowerCase()).filter((r) => r !== '')];
  return {
    name, label, tags, roles,
    description,
    signature: signatureOf(f),
    search: `${name} ${label} ${description}`.toLowerCase(),
    data: f,
  };
}

/** The handler is resolved per function — `find({})` mixes Func, Script and DataQuery, so one
 * cached handler would stamp the first row's icon on every row. Handler output is foreign plugin
 * DOM: a throwing handler degrades to no icon; the base handler's caption text is not an icon. */
function funcIcon(f: DG.Func): HTMLElement | null {
  try {
    const handler = DG.ObjectHandler.forEntity(f);
    const icon = handler ? handler.renderIcon(f) : null;
    return icon !== null && (icon.textContent ?? '').trim() === '' ? icon : null;
  } catch (_) {
    return null;
  }
}

export function functionsBrowser(options: DgFunctionsBrowserOptions = {}): FunctionsBrowser {
  let funcs = DG.Func.find({});
  const filter = options.filter;
  if (filter)
    funcs = funcs.filter(filter);
  const ignore = options.ignorePackages ?? [];
  if (ignore.length > 0)
    funcs = funcs.filter((f) => !ignore.some((p) => f.nqName.startsWith(`${p}:`)));
  if (options.scalarOnly)
    funcs = funcs.filter(hasScalarOutput);
  const items = funcs.map(funcItem);
  // ignorePunctuation: bare localeCompare fronts the registry's `_helper` / `/molecule.json`
  // registrations; the Dart view's ordering ignores the punctuation too
  const collator = new Intl.Collator('en', {ignorePunctuation: true, numeric: true});
  items.sort((a, b) => collator.compare(a.label, b.label));

  const onChanged = options.onChanged;
  return new FunctionsBrowser({
    ...options,
    items,
    roles: DG.functionRoles.map((r) => ({role: r.role, description: r.description})),
    runAction: (item) => (item.data as DG.Func).prepare().edit(),
    icon: (item) => funcIcon(item.data as DG.Func),
    onChanged: (item) => {
      if (item !== null && options.setCurrentObject !== false)
        grok.shell.o = item.data;
      onChanged?.(item);
    },
  });
}
