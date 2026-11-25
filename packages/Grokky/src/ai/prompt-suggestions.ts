import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as jsyaml from 'js-yaml';

import {_package} from '../package';

type Scope = 'column' | 'view' | 'global';
type Context = 'powerSearch' | 'panel';
export type SearchProvider = 'Help' | 'Execute' | 'Query';

export interface Suggestion {
  label?: string;
  /** The prompt sent to the assistant. Omitted for `action` suggestions handled client-side. */
  prompt?: string;
  /** A client-side action key (see [SUGGESTION_ACTIONS]) — runs immediately instead of prompting the AI. */
  action?: string;
  prefer?: Record<string, string>;
  /** Packages whose functions must be registered for this prompt to work (e.g. [Chem]). */
  requires?: string[];
}

export interface Block {
  context: Context;
  scope: Scope;
  key: string;
  match?: {semType?: string; type?: string};
  label: string;
  icon?: string;
  color?: string;
  requires?: string[];
  provider?: SearchProvider;
  suggestions: Suggestion[];
}

// TODO: LLM-curated suggestions vs current human + view-agnostic heuristics?
const SLOT = /\{(\w+)\.(\w+)}/g;
let _cache: Block[] | null = null;
const MAX_GROUPING_CATEGORIES = 50;

enum SlotFilter {
  Numerical = 'isNumerical',
  Categorical = 'isCategorical',
  DateTime = 'isDateTime',
  SemType = 'semType',
}

type SlotColumns = (filter: SlotFilter, semType?: string) => DG.Column[];

const SLOT_COLUMNS: Record<SlotFilter, (cols: DG.ColumnList, semType?: string) => Iterable<DG.Column>> = {
  [SlotFilter.Numerical]: (cols) => cols.numericalNoDateTime,
  [SlotFilter.Categorical]: (cols) => cols.categorical,
  [SlotFilter.DateTime]: (cols) => cols.dateTime,
  [SlotFilter.SemType]: (cols, semType) => semType ? cols.bySemTypeAll(semType) : [],
};

const isSlotFilter = (filter: string): filter is SlotFilter => filter in SLOT_COLUMNS;

function isGrouping(c: DG.Column): boolean {
  const n = c.categories.length;
  return n > 1 && n <= MAX_GROUPING_CATEGORIES;
}

function slotColumns(df: DG.DataFrame): SlotColumns {
  const lists = new Map<string, DG.Column[]>();
  return (filter, semType) => {
    const key = filter === SlotFilter.SemType ? `${filter}:${semType}` : filter;
    let list = lists.get(key);
    if (!list) {
      const first: DG.Column[] = [];
      const rest: DG.Column[] = [];
      for (const c of SLOT_COLUMNS[filter](df.columns, semType))
        (filter === SlotFilter.Categorical && !isGrouping(c) ? rest : first).push(c);
      lists.set(key, list = first.concat(rest));
    }
    return list;
  };
}

const _pkgAvailable = new Set<string>();

function requiresMet(req?: string[]): boolean {
  return !req || req.every((p) => {
    if (!_pkgAvailable.has(p) && DG.Func.find({package: p}).length > 0)
      _pkgAvailable.add(p);
    return _pkgAvailable.has(p);
  });
}

async function load(): Promise<Block[]> {
  if (!_cache)
    _cache = jsyaml.load(await _package.files.readAsText('suggestions.yaml')) as Block[];
  return _cache;
}

function resolveSuggestion(s: Suggestion, columns: SlotColumns | null, semType?: string): Suggestion | null {
  const pickedCols: Record<string, DG.Column> = {};
  const taken = new Set<DG.Column>();

  const slotSource = `${s.label ?? ''}\n${s.prompt ?? ''}`;
  for (const [, name, filter] of slotSource.matchAll(SLOT)) {
    if (name in pickedCols) continue;
    if (!columns || !isSlotFilter(filter)) return null;
    const candidates = columns(filter, semType).filter((c) => !taken.has(c));
    if (!candidates.length) return null;

    const prefer = s.prefer?.[name];
    const re = prefer ? new RegExp(prefer.match(/name~="([^"]+)"/)?.[1] ?? prefer, 'i') : null;
    const pick = (re && candidates.find((c) => re.test(c.name))) || candidates[0];
    pickedCols[name] = pick;
    taken.add(pick);
  }

  const fill = (t: string): string => t
    .replace(/\{(\w+)((?:\.\w+)+)}/g, (_, n, path) => {
      if (!path.includes('.', 1))
        return pickedCols[n]?.name ?? '';
      let val: unknown = pickedCols[n];
      for (const key of path.slice(1).split('.'))
        val = (val as Record<string, unknown> | undefined)?.[key];
      return `${val ?? ''}`;
    });
  return {
    label: s.label ? fill(s.label) : undefined,
    prompt: s.prompt ? fill(s.prompt) : undefined,
    action: s.action,
  };
}

export async function resolveScopes(context: Context, view: DG.ViewBase | null = null): Promise<Block[]> {
  const df = view instanceof DG.TableView ? view.dataFrame ?? null : null;
  const columns = df ? slotColumns(df) : null;
  const out: Block[] = [];
  for (const b of await load()) {
    if (b.context !== context) continue;
    const applies = b.scope === 'global' ||
      (b.scope === 'view' && view?.type === b.match?.type) ||
      (b.scope === 'column' && !!df?.columns.bySemType(b.match?.semType ?? ''));
    if (!applies || !requiresMet(b.requires)) continue;

    const suggestions = b.suggestions
      .map((s) => requiresMet(s.requires) ? resolveSuggestion(s, columns, b.match?.semType) : null)
      .filter((s): s is Suggestion => !!s);
    if (suggestions.length)
      out.push({...b, suggestions});
  }
  return out;
}

export function showSuggestionsMenu(
  scopes: Block[],
  onSelect: (s: Suggestion, provider?: SearchProvider) => void,
  causedBy?: MouseEvent,
): void {
  const menu = DG.Menu.popup();
  for (const s of scopes) {
    const group = menu.group(s.label);
    for (const sg of s.suggestions)
      group.item(sg.label ?? sg.prompt ?? '', () => onSelect(sg, s.provider));
  }
  menu.show(causedBy ? {causedBy} : undefined);
}

// ---- Client-side actions -------------------------------------------------------------
// Some suggestions are ambiguous for the AI (e.g. "open a demo dataset"): resolving them
// with a quick in-panel choice is faster and more predictable than a round-trip to the assistant.

/** One selectable option in an in-panel choice block. */
export interface ChoiceOption {
  label: string;
  onSelect: () => void;
}

/** The panel surface an action renders into — it draws the choices as an assistant-style response. */
export interface SuggestionActionHost {
  addChoice(prompt: string | null, options: ChoiceOption[]): void;
  addNote(markdown: string): void;
}

const DEMO_DATASETS: {label: string; open: () => Promise<DG.DataFrame>}[] = [
  {label: 'Demographics (clinical)', open: () => grok.dapi.files.readCsv('System:DemoFiles/demog.csv')},
  {label: 'SPGI (small molecules)', open: () => grok.dapi.files.readCsv('System:DemoFiles/chem/SPGI.csv')},
  {label: 'Dose-response curves', open: () => grok.dapi.files.readCsv('System:DemoFiles/curves.csv')},
  {label: 'Beer reviews', open: () => grok.dapi.files.readCsv('System:DemoFiles/beer.csv')},
  {label: 'Random walk', open: async () => grok.data.demo.randomWalk()},
];

function openDemoDataset(d: {label: string; open: () => Promise<DG.DataFrame>}, host: SuggestionActionHost): void {
  d.open()
    .then((df) => {
      grok.shell.addTableView(df);
      host.addNote(`Opened the **${d.label}** table.`);
    })
    .catch((e) => grok.shell.error(`Couldn't open ${d.label}: ${e}`));
}

const SUGGESTION_ACTIONS: Record<string, (host: SuggestionActionHost) => void> = {
  openDemoDataset: (host) => host.addChoice(
    'Which demo dataset would you like to open?',
    DEMO_DATASETS.map((d) => ({label: d.label, onSelect: () => openDemoDataset(d, host)})),
  ),
};

/** Runs a suggestion's client-side action into the given panel host, if it has one. Returns whether it was handled. */
export function runSuggestionAction(action: string | undefined, host: SuggestionActionHost): boolean {
  const handler = action ? SUGGESTION_ACTIONS[action] : undefined;
  if (handler)
    handler(host);
  return !!handler;
}
