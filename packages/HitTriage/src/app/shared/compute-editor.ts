/* ComputeEditor: what a template calculates for every molecule, the descriptors and the compute functions, as one
   master-detail box: a row per function on the left (its checkbox, its name, how many descriptors are picked), the
   current row's settings on the right: the descriptors as a searchable tree whose groups say how many of theirs are
   picked, or a function's own form and "Re-run when campaign opens". "Reorder…" orders and hides the functions, kept
   where the compute dialog keeps them. The value is the compute config the compute dialog writes. Nothing here is
   HitTriage's but the box's size (`.hit-triage-compute-box` in css/hit-triage.css, which moves with it). For Chem MPO,
   which asks for one function and no descriptors, it still lacks the dialog's single choice (`singleSelect`). */
import * as DG from 'datagrok-api/dg';
import {BoolInput, ChecklistInput, Control, Dialog, MasterDetail, TextInput, VirtualTree, component, computed, div,
  divH, divV, link, signal, statusLine, untracked} from '@datagrok-libraries/u2';
import type {Signal, TreeNode} from '@datagrok-libraries/u2';
import {funcForm} from '@datagrok-libraries/u2/dg';
import type {ComputeFunctions, IComputeDialogResult, IDescriptorTree, IFunctionArgs, TemplateCompute}
  from '@datagrok-libraries/statistics/src/compute-functions/types';
import {HTQueryPrefix, HTScriptPrefix} from '@datagrok-libraries/statistics/src/compute-functions/consts';
import {getFuncPackageNameSafe, getSavedFunctionOrdering, setSavedFunctionOrdering}
  from '@datagrok-libraries/statistics/src/compute-functions/utils';
import '../../../css/hit-triage.css';

export interface ComputeEditorOptions {
  /** The functions a template may calculate (`discoverComputeFunctions`). */
  functions: ComputeFunctions;
  value: TemplateCompute;
  /** The descriptors to pick from; without it there is no Descriptors row and the value's descriptors stay as given. */
  descriptors?: () => Promise<IDescriptorTree>;
  /** Where the order of the functions and the hidden ones are kept, as {order, hidden} of their names. */
  orderKey: string;
  /** A checkbox per function that flags it to be calculated again whenever a campaign opens. */
  rerunOnOpen?: {label: string, tooltip: string};
}

const DESCRIPTORS = 'Descriptors';

interface Entry {
  /** `Descriptors`, `<package>:<name>`, `HTScript:<name>:<id>` or `HTQuery:<name>:<id>`: the compute dialog's keys. */
  key: string;
  /** The friendly name: what the row shows and what the order is kept by. */
  name: string;
  checked: Signal<boolean>;
  rerun: Signal<boolean>;
  /** The call holding the function's arguments; none for the descriptors. */
  call?: DG.FuncCall;
  /** The inputs a template does not keep: a function's table and column, a query's molecules. */
  skip: number;
}

type Kept = {args?: IFunctionArgs, rerunOnOpen?: boolean} | undefined;

export interface ComputeEditor {
  /** The descriptors are listed (or could not be): until then the ones picked read as given. */
  readonly ready: Promise<void>;
  /** A pick in the descriptor tree changed what the editor holds (the inputs tell their own edits). */
  readonly edited: Signal<boolean>;
  /** The compute config the editor holds, in the compute dialog's shape. */
  read(): TemplateCompute;
}

export function computeEditor(options: ComputeEditorOptions): Control & ComputeEditor {
  return component((self) => {
    const edited = signal(false);
    const tree = signal<IDescriptorTree | null>(null);
    const treeError = signal<string | null>(null);
    const given = options.value.descriptors;
    const picked = signal(new Set(given?.args ?? []));
    const descriptorEntry: Entry | undefined = options.descriptors === undefined ? undefined : {key: DESCRIPTORS,
      name: DESCRIPTORS, skip: 0, rerun: signal(false),
      checked: signal(!!given?.enabled && (given?.args?.length ?? 0) > 0)};
    const entries = [...descriptorEntry ? [descriptorEntry] : [], ...computeEntries(options.functions, options.value)];
    const order = signal(getSavedFunctionOrdering(options.orderKey));
    const ready = options.descriptors === undefined ? Promise.resolve() : options.descriptors().then((t) => {
      tree.value = t;
    }, (e) => {
      treeError.value = e instanceof Error ? e.message : String(e);
    });

    // built now, not when first shown: a function's form writes its computed defaults into the call read() keeps
    const details = new Map<string, Control | HTMLElement>();
    for (const e of entries)
      details.set(e.key, e.call === undefined ? descriptorsDetail() : functionDetail(e));
    const box = new MasterDetail<Entry>({items: computed(() => shown()), keyOf: (e) => e.key, label: 'Functions',
      row: (e) => ({label: e.name, check: e.checked, checkTooltip: 'Toggle calculation of this function',
        hint: e.call !== undefined ? undefined : computed(() => {
          const n = descriptorArgs().length;
          return e.checked.value && n > 0 ? `${n} selected` : '';
        })}),
      detail: (e) => details.get(e.key)!,
      disabled: (e) => e.call === undefined ? null :
        computed(() => e.checked.value ? null : `Check ${e.name} to change its settings`)});
    box.root.classList.add('hit-triage-compute-box');
    const root = divV([divH([link('Reorder…', () => reorder())]), box], {gap: 's'});
    return {root, ready, edited, read};

    function read(): TemplateCompute {
      const hidden = new Set(order.peek().hidden);
      const result: IComputeDialogResult = {descriptors: descriptorEntry || !given?.enabled ? [] : given.args ?? [],
        externals: {}, scripts: {}, queries: {}, ...options.rerunOnOpen ? {rerunOnOpen: {}} : {}};
      for (const e of entries) {
        if (!e.checked.peek() || hidden.has(e.name))
          continue;
        if (e.call === undefined) {
          result.descriptors = untracked(() => descriptorArgs());
          continue;
        }
        const args: IFunctionArgs = {};
        for (const [name, value] of Object.entries(e.call.inputs).slice(e.skip))
          args[name] = value;
        const target = e.key.startsWith(`${HTScriptPrefix}:`) ? result.scripts! :
          e.key.startsWith(`${HTQueryPrefix}:`) ? result.queries! : result.externals;
        target[e.key] = args;
        if (result.rerunOnOpen && e.rerun.peek())
          result.rerunOnOpen[e.key] = true;
      }
      return templateCompute(result);
    }

    /** The entries in the saved order of their names, the ones it does not name after them, the hidden ones out. */
    function shown(): Entry[] {
      const {order: names, hidden} = order.value;
      const rank = (e: Entry) => names.includes(e.name) ? names.indexOf(e.name) : names.length;
      return entries.map((e, at) => ({e, at})).filter(({e}) => !hidden.includes(e.name))
        .sort((a, b) => rank(a.e) - rank(b.e) || a.at - b.at).map(({e}) => e);
    }

    function descriptorsDetail(): HTMLElement {
      const search = new TextInput({inline: true, search: true, name: 'descriptor-search',
        placeholder: 'Search descriptors'});
      const view = new VirtualTree<string>({itemHeight: 24, onCheck: (node, checked) => check(node, checked)});
      const status = div([statusLine(computed(() => treeError.value === null ? null :
        `Could not list the descriptors: ${treeError.value}`), 'Listing the descriptors…')]);
      self.effect(() => {
        const query = search.value.trim().toLowerCase();
        const roots = nodes(query);
        view.setRoots(roots);
        if (query !== '')
          view.expanded.value = new Set(roots.map((r) => r.id));
      });
      self.effect(() => {
        const listed = tree.value !== null;
        view.root.hidden = !listed;
        status.hidden = listed;
      });
      return divV([search.root, status, view.root], {gap: 's', cls: 'hit-triage-compute-descriptors'});
    }

    /** The function's own form, its table and column hidden, and the re-run flag last on the same label column. */
    function functionDetail(e: Entry): Control {
      const hide = Object.fromEntries(e.call!.func.inputs.slice(0, e.skip).map((p) => [p.name, {visible: false}]));
      const form = funcForm(e.call!,
        {showTableSelectors: false, skipTableAutoFill: true, editors: 'auto', overrides: hide});
      const rerun = options.rerunOnOpen;
      if (rerun) {
        form.add(new BoolInput({label: rerun.label, name: `rerun-${e.key}`, bind: e.rerun,
          tooltipText: rerun.tooltip}));
      }
      return form;
    }

    /** The descriptor groups, each with the descriptors that match `query` and how many of all its own are picked. */
    function nodes(query: string): TreeNode<string>[] {
      const groups = tree.value;
      const names = picked.value;
      if (groups === null)
        return [];
      const roots: TreeNode<string>[] = [];
      for (const [group, {descriptors}] of Object.entries(groups)) {
        const matching = query === '' ? descriptors : descriptors.filter((d) => d.name.toLowerCase().includes(query));
        if (matching.length === 0)
          continue;
        const on = descriptors.filter((d) => names.has(d.name)).length;
        roots.push({id: `group:${group}`, label: `${group} · ${on} of ${descriptors.length}`, checked: on > 0,
          children: matching.map((d) => ({id: `${group}/${d.name}`, label: d.name, tooltip: d.description,
            data: d.name, checked: names.has(d.name)}))});
      }
      return roots;
    }

    /** A descriptor, or a group's descriptors the search shows, picked or dropped; picking turns the descriptors on. */
    function check(node: TreeNode<string>, checked: boolean): void {
      const names = node.data !== undefined ? [node.data] : (node.children as TreeNode<string>[]).map((c) => c.data!);
      const next = new Set(picked.peek());
      for (const name of names) {
        if (checked)
          next.add(name);
        else
          next.delete(name);
      }
      picked.value = next;
      edited.value = true;
      if (checked)
        descriptorEntry!.checked.value = true;
    }

    /** The picked descriptors in the tree's order; as given while the tree is not listed. */
    function descriptorArgs(): string[] {
      const groups = tree.value;
      const names = picked.value;
      if (groups === null)
        return [...names];
      return Object.values(groups).flatMap((g) => g.descriptors.map((d) => d.name)).filter((n) => names.has(n));
    }

    function reorder(): void {
      self.build(() => {
        const names = ordered(entries.map((e) => e.name), order.peek().order);
        const hidden = order.peek().hidden;
        const functions = new ChecklistInput({name: 'functions', orderable: true, items: names,
          value: names.filter((n) => !hidden.includes(n))});
        new Dialog('Order and hide functions')
          .add(functions)
          .onOK(() => {
            const kept = functions.peek();
            order.value = {order: kept, hidden: names.filter((n) => !kept.includes(n))};
            setSavedFunctionOrdering(order.peek(), options.orderKey);
          })
          .show();
      });
    }
  }, {name: 'compute-editor'});
}

/** The functions, scripts and queries a template may calculate, each with what `value` keeps of it; one that cannot
 * be prepared is left out, as the compute dialog leaves it. */
function computeEntries(functions: ComputeFunctions, value: TemplateCompute): Entry[] {
  const entries: Entry[] = [];
  const add = (f: DG.Func, key: () => string, kept: () => Kept, skip: number) => {
    try {
      const k = kept();
      entries.push({key: key(), name: f.friendlyName ?? f.name, skip, checked: signal(k !== undefined),
        rerun: signal(!!k?.rerunOnOpen), call: f.prepare(k?.args ?? {})});
    }
    catch (e) {
      console.error(e);
    }
  };
  const onColumn = (f: DG.Func) => f.inputs.length >= 2 && f.inputs[0].propertyType === 'dataframe' &&
    f.inputs[1].propertyType === 'column';
  for (const f of functions.functions.filter(onColumn)) {
    add(f, () => `${f.package.name}:${f.name}`, () => value.functions?.find((t) => t.name === f.name &&
      t.package === getFuncPackageNameSafe(f)), 2);
  }
  for (const s of functions.scripts.filter(onColumn))
    add(s, () => `${HTScriptPrefix}:${s.name ?? ''}:${s.id}`, () => value.scripts?.find((t) => t.id === s.id), 2);
  for (const q of functions.queries.filter((q) => q.inputs.length > 0 && q.inputs[0].propertyType === 'list'))
    add(q, () => `${HTQueryPrefix}:${q.name ?? ''}:${q.id}`, () => value.queries?.find((t) => t.id === q.id), 1);
  return entries;
}

/** `names`, each once, in the saved order, the ones it does not name after them. */
function ordered(names: string[], order: string[]): string[] {
  const unique = [...new Set(names)];
  const known = order.filter((n) => unique.includes(n));
  return [...known, ...unique.filter((n) => !known.includes(n))];
}

/** The compute config a template keeps, from the compute dialog's answer; `rerunOnOpen` is written only where the
 * dialog offered it. */
export function templateCompute(result: IComputeDialogResult): TemplateCompute {
  const rerun = (key: string) => result.rerunOnOpen ? {rerunOnOpen: !!result.rerunOnOpen[key]} : {};
  const byId = (entries: {[key: string]: IFunctionArgs} | undefined, prefix: string) =>
    Object.entries(entries ?? {}).filter(([key]) => key.startsWith(prefix) && key.split(':').length === 3)
      .map(([key, args]) => {
        const [, name, id] = key.split(':');
        return {name, id, args, ...rerun(key)};
      });
  return {
    descriptors: {enabled: result.descriptors.length > 0, args: result.descriptors},
    functions: Object.entries(result.externals).map(([key, args]) => {
      const [pkg, name] = key.split(':');
      return {name, package: pkg, args, ...rerun(key)};
    }),
    scripts: byId(result.scripts, HTScriptPrefix),
    queries: byId(result.queries, HTQueryPrefix),
  };
}
