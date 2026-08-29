/* One `u2-viewer-<kebab-name>` tag per platform viewer descriptor (V3, VP-9), registered from the
   live list by an explicit call — never at import: the gallery, the manifest script and the
   headless suite import the core registry, which must stay platform-free. */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {Scope} from '../../core/scope.js';
import {computed, Signal} from '../../core/signals.js';
import {Component} from '../../core/component.js';
import type {LiveOption} from '../../core/input-base.js';
import {Dialog} from '../../components/containers/dialog.js';
import {ChoiceInput} from '../../components/inputs/choice-input.js';
import {Registry, registry as globalRegistry} from '../../spec/registry.js';
import type {ComponentMeta, DesignerAction, SpecPropMeta} from '../../spec/registry.js';
import {registerAll} from '../../spec/registrations.js';
import {specWarn} from '../../spec/spec.js';
import type {SpecNode} from '../../spec/spec.js';
import {columnInput} from '../inputs/pickers.js';
import {functionsBrowser} from '../entities/functions-browser.js';
import {FunctionInput} from '../inputs/function-input.js';
import {viewerControl} from './viewer-control.js';

const api = globalThis as {grok_Property_Get_PropertySubType?: (dart: unknown) => string | null};

export const kebab = (name: string): string =>
  name.toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-|-$/g, '');

const TABLE: SpecPropMeta = {name: 'table', type: 'dataframe', bindable: true,
  description: 'The DataFrame this viewer shows — bind it to a source (`$.orders`).'};
const LOOK: SpecPropMeta = {name: 'look', type: 'object', category: 'Misc',
  description: 'Any further look keys (the I*Settings shape); a named prop wins over the same key here.'};

const COPIED = ['friendlyName', 'description', 'category', 'semType', 'min', 'max', 'step',
  'inputType', 'editor', 'format', 'units', 'showSlider', 'showPlusMinus'] as const;

/** Dart-internal structures the platform dumps as JSON, never a look a user edits. */
const INTERNAL_TYPES = new Set(['lookandfeel', 'histogramlook', 'gridlook', 'gridcellstyle']);
/** Platform script hooks (`onClick`, `initializationFunction`, `onInitializedScript`): not look. */
const SCRIPT_HOOK = /^on[A-Z]|Function$|Script$/;

/** `propertyType` is deliberately not copied: `checkProp` reads `propertyType ?? type`, and the
 * platform's `num`/`list` would fail its type check. Nor is `nullable`: the platform's panel never
 * marks a look prop required, and a copied `false` would turn every untouched empty prop red.
 * Read field by field — never spread a `DG.Property`, every field of which is a prototype getter. */
export function toSpecProp(p: DG.Property): SpecPropMeta {
  const out: SpecPropMeta = {name: p.name, type: specType(p), bindable: true, twoWay: p.set != null};
  for (const key of COPIED) {
    const value = p[key];
    if (value != null)
      (out as unknown as Record<string, unknown>)[key] = value;
  }
  const choices = p.choices;
  if (choices != null && choices.length > 0)
    out.choices = choices;
  return out;
}

function specType(p: DG.Property): string {
  const type = p.propertyType;
  switch (type) {
    case 'string': case 'int': case 'double': case 'bool': return type;
    case 'num': return 'double';
    case 'list': return subTypeOf(p) === 'string' ? 'string_list' : 'object';
    default: return 'object';
  }
}

const subTypeOf = (p: DG.Property): string | null => api.grok_Property_Get_PropertySubType?.(p.dart) ?? null;

/** The descriptor's user look: not its own `table`/`viewer` (the u2 `table` bind covers them), not
 * the `Events` category (platform script hooks, not u2 events), not a Dart-internal structure. */
const isUserLook = (p: DG.Property): boolean =>
  p.userEditable !== false && p.name !== 'table' && p.name !== 'viewer' && p.category !== 'Events' &&
  !INTERNAL_TYPES.has(p.propertyType) && !(p.propertyType === 'string' && SCRIPT_HOOK.test(p.name));

export const VIEWER_USAGE: Record<string, string> = {
  'Grid': 'Tabular detail over a frame — the master list of a master-detail. Bind `table` to a source; ' +
    'rows, selection and the current row sync through the frame, so a form bound to ' +
    '`$.orders.currentRow.*` follows a click with no wiring. Not for a list of arbitrary items — ' +
    'that is `VirtualList`/`VirtualGrid`.',
  'Filters': 'The filter panel over a frame; filtering reaches every other viewer on the same table ' +
    'through the frame. Without `filters` it shows every column; set `filters` (the FilterState list) ' +
    'or use *Add filter for column…* to choose.',
  'Form': 'Field-per-column view of the current row, in the viewer\'s own layout. Prefer `u2-form` with ' +
    'bound inputs when you want u2 inputs and validation; prefer this for a quick record view with the ' +
    'platform\'s field editors.',
  'Scatter plot': 'Two numeric columns, optional color/size/marker columns. Bind ' +
    '`xColumnName`/`yColumnName` to choice inputs over `$.orders.columns` for a configurable plot.',
  'Bar chart': 'Category counts or an aggregate per category (`splitColumnName`, `valueColumnName`, ' +
    '`valueAggrType`); horizontal bars by default.',
  'Histogram': 'Distribution of one numeric column (`valueColumnName`); doubles as a range filter when ' +
    '`filteringEnabled`.',
  'Line chart': 'One or more numeric series (`yColumnNames`) over an x column, typically time.',
  'Pie chart': 'Share per category (`categoryColumnName`); for more than ~7 categories a bar chart reads better.',
  'Tile Viewer': 'Cards, one per row, laid out in a grid — the browse view over a frame; for cards over ' +
    'an item array use u2\'s own `VirtualGrid`.',
};

/** 'Add filter for column…' on `u2-viewer-filters` (VP-11): the column picked in a dialog joins the
 * `filters` prop as one more `FilterState`; the set-prop's rerender is what shows it (P14). The
 * list is seeded from the live group — the platform turns `columnNames` into `filters` at build,
 * so the document's `filters` alone would drop the panes it already shows. Read through the
 * serialized look: `props.filters` answers Dart maps, which are not JSON. */
function filterSeed(node: SpecNode, built: DG.Viewer): DG.FilterState[] {
  const live = built.getOptions().look.filters as DG.FilterState[] | undefined;
  return [...(live != null && live.length > 0 ? live : node.props?.filters as DG.FilterState[] | undefined ?? [])];
}

/** One choice, OK'd in a modal dialog: the pick, or null on Cancel — and when `owner` (the viewer's
 * scope) goes away first, which cancels the dialog with the node or its view. */
function ask(title: string, owner: Scope | undefined, make: () => ChoiceInput): Promise<string | null> {
  return new Promise((resolve) => {
    const scope = new Scope();
    const finish = (picked: string | null) => {
      scope.dispose();
      resolve(picked);
    };
    const cancel = () => finish(null);
    owner?.own(cancel);
    scope.own(() => owner?.disown(cancel));
    try {
      Scope.runWith(scope, () => {
        const input = make();
        new Dialog(title)
          .add(input)
          .onOK(() => finish(input.value.peek()))
          .onCancel(cancel)
          .okEnabled(computed(() => input.value.value !== null))
          .show({modal: true, width: 360});
      });
    } catch (e) {
      scope.dispose();
      throw e;
    }
  });
}

export const ADD_FILTER: DesignerAction = {
  name: 'Add filter for column…',
  icon: 'filter',
  produce: async (node, built) => {
    const viewer = built instanceof DG.Viewer ? built : null;
    if (viewer?.dataFrame == null || viewer.dataFrame.columns.names().length === 0) {
      grok.shell.warning(`${ADD_FILTER.name}: no table to pick a column from — bind \`table\` to a source first`);
      return null;
    }
    const df = viewer.dataFrame;
    const filters = filterSeed(node, viewer);
    const filtered = new Set(filters.map((f) => f.column));
    const name = await ask('Add filter', viewer.scope,
      () => columnInput('Column', df, {filter: (c) => !filtered.has(c.name)}));
    const col = name === null ? null : df.columns.byName(name);
    return col === null ? null : {op: 'set-prop', name: 'filters',
      value: [...filters, {type: col.isNumerical ? 'histogram' : 'categorical', column: name}]};
  },
};

/** 'Remove filter…' on `u2-viewer-filters`: one of the DOCUMENT's `filters` entries picked in a
 * dialog leaves the list — the key itself when it was the last; the panes the platform built from
 * `columnNames` are not the document's to remove. */
export const REMOVE_FILTER: DesignerAction = {
  name: 'Remove filter…',
  icon: 'filter',
  produce: async (node, built) => {
    const filters = node.props?.filters as DG.FilterState[] | undefined ?? [];
    if (filters.length === 0) {
      grok.shell.warning(`${REMOVE_FILTER.name}: the form names no filters — "${ADD_FILTER.name}" adds one`);
      return null;
    }
    const picked = await ask('Remove filter', Component.is(built) ? built.scope : undefined,
      () => new ChoiceInput({label: 'Filter', nullable: true,
        items: filters.map((f, i) => ({value: String(i), label: f.column ?? `#${i + 1}`}))}));
    if (picked === null)
      return null;
    const rest = filters.filter((_, i) => i !== Number(picked));
    return {op: 'set-prop', name: 'filters', value: rest.length > 0 ? rest : undefined};
  },
};

function viewerMeta(d: DG.WidgetDescriptor, tag: string): ComponentMeta {
  return {
    tag,
    label: d.name,
    icon: () => d.createIcon() as HTMLElement,
    category: 'Viewers',
    appearance: false,
    description: d.description,
    usage: VIEWER_USAGE[d.name],
    props: [TABLE, LOOK, ...d.properties.filter(isUserLook).map(toSpecProp)],
    events: [...new Set(['d4-data-event', ...d.events.map((e) => e.name)])],
    create: (props) => viewerControl(d.name, props),
    example: {tag, bind: {table: '$.orders'}},
    designerActions: d.name === 'Filters' ? [ADD_FILTER, REMOVE_FILTER] : undefined,
  };
}

const registered = new WeakSet<Registry>();

/** Every descriptor of `DG.WidgetDescriptor.getDescriptors()` as a `u2-viewer-*` tag; once per
 * registry; a package-internal one (`_makeInspectorPanel`) is passed over, and a tag already
 * taken — by a collision of kebab names, or an earlier registration — is skipped with one warning.
 * Call it from an app or autostart function, never at module scope: the descriptors exist only
 * once the platform is up. */
export function registerPlatformViewers(reg: Registry = globalRegistry): void {
  if (registered.has(reg))
    return;
  const descriptors = DG.WidgetDescriptor.getDescriptors();
  registered.add(reg);
  const taken = new Map<string, string>();
  for (const d of descriptors) {
    if (d.name.startsWith('_'))
      continue;
    const tag = `u2-viewer-${kebab(d.name)}`;
    if (reg.get(tag) !== undefined) {
      const by = taken.get(tag);
      specWarn(by === d.name ? `${tag}: "${d.name}" is already registered` :
        `${tag}: "${d.name}" skipped — the tag is taken by "${by ?? 'an earlier registration'}"`);
      continue;
    }
    reg.register(viewerMeta(d, tag));
    taken.set(tag, d.name);
  }
}

const FB_TAG = 'u2-functions-browser';

function lit<T>(x: unknown): T | undefined {
  return x instanceof Signal ? (x as Signal<T>).peek() : x as T | undefined;
}

/** The function registry browser (see the core `FunctionsBrowser`). Platform-side, like the
 * viewer tags — it needs `DG.Func.find`, so it stays out of the platform-free manifest. The
 * declarative subset only; the function-typed options stay a code-level escape hatch. */
function functionsBrowserMeta(): ComponentMeta {
  return {
    tag: FB_TAG,
    label: 'Functions browser',
    category: 'Data',
    description: 'Searchable, tag/role-filterable virtualized list of every registered function.',
    usage: 'The function registry browser: plain terms filter by name, `#tag` / `@role` terms ' +
      'check the filter panes and vice versa. Bind `search` two-way to drive it from outside; ' +
      '`change` fires on selection (which also becomes the shell\'s current object), `activate` ' +
      'on double-click or Enter.',
    create: (props) => functionsBrowser({
      search: lit<string>(props.search),
      tags: lit<string[]>(props.tags),
      visibleTags: lit<string[]>(props.visibleTags),
      ignoreTags: lit<string[]>(props.ignoreTags),
      ignorePackages: lit<string[]>(props.ignorePackages),
      scalarOnly: lit<boolean>(props.scalarOnly),
      showSearch: props.showSearch as LiveOption<boolean> | undefined,
      showTags: props.showTags as LiveOption<boolean> | undefined,
      showSignature: props.showSignature as LiveOption<boolean> | undefined,
      showRunButton: props.showRunButton as LiveOption<boolean> | undefined,
    }),
    props: [
      {name: 'search', type: 'string', bindable: true, twoWay: true,
        description: 'The query: plain terms, `#tag` and `@role`.'},
      {name: 'showSearch', type: 'bool', bindable: true, description: 'Shows the search bar.'},
      {name: 'showTags', type: 'bool', bindable: true,
        description: 'Shows the roles/tags filter panes.'},
      {name: 'showSignature', type: 'bool', bindable: true,
        description: 'Shows the `(inputs) : outType` signature on every row.'},
      {name: 'showRunButton', type: 'bool', bindable: true,
        description: 'Shows the per-row play icon.'},
      {name: 'tags', type: 'string_list', description: 'Pre-checked tags.'},
      {name: 'visibleTags', type: 'string_list',
        description: 'Restricts the tags pane to these tags.'},
      {name: 'ignoreTags', type: 'string_list',
        description: 'Tags left out of the tags pane; res, converters, internal by default.'},
      {name: 'ignorePackages', type: 'string_list',
        description: 'Packages whose functions are left out.'},
      {name: 'scalarOnly', type: 'bool',
        description: 'Only functions with one scalar (or dynamic) output.'},
    ],
    events: ['change', 'activate'],
    example: {tag: FB_TAG, name: 'fb', props: {showSignature: false}},
  };
}

const FI_TAG = 'u2-function-input';

/** The FunctionName picker (see `FunctionInput`). Platform-side like the browser it opens —
 * it needs `DG.Func.find`. The declarative subset only; `filter` stays a code-level option. */
function functionInputMeta(): ComponentMeta {
  return {
    tag: FI_TAG,
    label: 'Function input',
    category: 'Inputs',
    description: 'Namespace-qualified function name, picked from the FunctionsBrowser in a popup.',
    usage: 'The value is the nqName string (`Chem:SmilesToMw`) — the `FunctionName` semantic ' +
      'type. A row click, Enter or double-click in the popup commits and closes.',
    create: (props) => {
      const value = props.value;
      const bound = value instanceof Signal;
      return new FunctionInput({
        label: props.label as LiveOption<string> | undefined,
        name: lit<string>(props.name),
        nullable: lit<boolean>(props.nullable),
        placeholder: lit<string>(props.placeholder),
        scalarOnly: lit<boolean>(props.scalarOnly),
        ignorePackages: lit<string[]>(props.ignorePackages),
        value: bound ? undefined : value as string | undefined,
        bind: bound ? value as Signal<string> : undefined,
      });
    },
    props: [
      {name: 'label', type: 'string', bindable: true},
      {name: 'name', type: 'string', description: 'Stable key for forms and dumps; defaults to the label.'},
      {name: 'value', type: 'string', bindable: true, twoWay: true,
        description: 'The namespace-qualified function name.'},
      {name: 'nullable', type: 'bool', description: 'Backspace/Delete clears the value.'},
      {name: 'placeholder', type: 'string', description: 'Shown while nothing is picked.'},
      {name: 'scalarOnly', type: 'bool',
        description: 'Only functions with one scalar (or dynamic) output.'},
      {name: 'ignorePackages', type: 'string_list',
        description: 'Packages whose functions are left out of the popup.'},
    ],
    events: ['input', 'change'],
    example: {tag: FI_TAG, props: {label: 'Function', value: 'Chem:SmilesToMw'}},
  };
}

/** The platform registry: every core tag plus every viewer tag plus the platform-side controls. */
export function registerPlatformComponents(reg: Registry = globalRegistry): void {
  registerAll(reg);
  registerPlatformViewers(reg);
  if (reg.get(FB_TAG) === undefined)
    reg.register(functionsBrowserMeta());
  if (reg.get(FI_TAG) === undefined)
    reg.register(functionInputMeta());
}
