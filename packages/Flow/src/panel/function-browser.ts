/* eslint-disable max-len */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {FuncInfo, getRegisteredFuncs, isWorkflowFunc, VIEWER_NODE_TYPES} from '../rete/node-factory';
import {tid, setTid} from '../utils/test-ids';
import {getFilesBrowser} from '../utils/files-browser-tree';
import {categorizeBySignature, domainCategory} from '../types/type-map';
import {getFavorites, isFavorite, toggleFavorite, onFavoritesChanged} from './favorites';
import {supportedUploadExtensions} from '../utils/uploaded-files';

/** Whether a function's output is a widget (→ the Widgets pane, not a category). */
export function funcOutputsWidget(f: FuncInfo): boolean {
  try {
    return f.func.outputs.some((o) => String(o.propertyType) === 'widget');
  } catch {
    return false;
  }
}

export type GroupByMode = 'category' | 'role' | 'tags' | 'package';

/** The task-oriented categories, in the order a scientist meets a pipeline:
 *  get data in → combine → reshape → derive columns → compute values →
 *  visualize → everything else. Used both to bucket and to order the tree. */
export const FUNC_CATEGORIES = [
  'Data Sources',
  'Workflows',
  'Combine Tables',
  'Transform Tables',
  'Column Operations',
  'Compute Values',
  'Visualize',
  'Cheminformatics',
  'Bioinformatics',
  'Other',
] as const;
export type FuncCategory = (typeof FUNC_CATEGORIES)[number];

/** Bucket a function by **what it does**, derived from its input/output
 *  signature (and role for viewers) — delegates to the shared
 *  `categorizeBySignature` (the same logic that colors node title bars). The
 *  key distinctions, validated against the live catalog:
 *  - **Data Sources** produce a table from *no* table input (OpenFile, DB
 *    queries, generators). A join is NOT a data source — it consumes tables.
 *  - **Combine Tables** take ≥2 tables (join / union / link / react).
 *  - **Transform Tables** take exactly one table and return a table or mutate
 *    it in place (filter, aggregate, pivot, sort, …).
 *  - **Column Operations** emit a column / column list (AddNewColumn, descriptors).
 *  - **Compute Values** emit a scalar (statistics, math, text, predicates).
 *  - **Visualize** emit a viewer / view / widget / graphics.
 *  - **Cheminformatics** / **Bioinformatics** — grouped by source package
 *    (domain wins over the signature category), so a scientist finds all
 *    chem/bio steps together regardless of what they do.
 *  - **Workflows** — saved flows (`DG.Script`, language `flow`); checked
 *    before everything else since their signature reads as a data source.
 *  - **Other** — expression/filter builders and dynamic helpers. */
export function categorizeFunc(func: DG.Func, role: string | null, packageName?: string): FuncCategory {
  if (isWorkflowFunc(func)) return 'Workflows';
  const inputTypes = func.inputs.map((i) => String(i.propertyType));
  const domain = domainCategory(packageName, inputTypes);
  if (domain) return domain;
  return categorizeBySignature(
    inputTypes,
    func.outputs.map((o) => String(o.propertyType)),
    role) as FuncCategory;
}

/** For each domain section, the flagship package whose functions lead the list
 *  (the rest of that domain's packages follow). */
export const DOMAIN_PRIMARY_PACKAGE: Record<string, string> = {
  Cheminformatics: 'Chem',
  Bioinformatics: 'Bio',
};

/** The "most-used" operations a scientist reaches for first inside a domain
 *  section — matched as lowercased substrings against the function name (lower
 *  index = higher up). Everything unlisted sorts after these, alphabetically —
 *  so descriptors / properties / similarity lead each domain list. */
export const DOMAIN_PRIORITY_KEYWORDS: Record<string, string[]> = {
  Cheminformatics: [
    'descriptor', 'propert', 'fingerprint', 'similar', 'substructure',
    'diversit', 'cluster', 'scaffold', 'r-group', 'rgroup', 'mcs',
    'activity', 'toxic', 'risk', 'elemental', 'mutate', 'curate', 'convert',
  ],
  Bioinformatics: [
    'descriptor', 'propert', 'composition', 'similar', 'sequence',
    'align', 'msa', 'logo', 'atomic', 'translate', 'helm', 'convert',
  ],
};

/** Order a domain section's items: flagship-package functions first (Chem under
 *  Cheminformatics, Bio under Bioinformatics), then the rest; within each group
 *  the curated "most-used" operations (descriptors, properties, …) lead, then
 *  alphabetical by package, then by display name. Pure — returns a new array. */
export function orderDomainSection(items: FuncInfo[], category: string): FuncInfo[] {
  const primary = DOMAIN_PRIMARY_PACKAGE[category];
  const keywords = DOMAIN_PRIORITY_KEYWORDS[category] ?? [];
  const priorityRank = (info: FuncInfo): number => {
    const hay = `${info.name} ${info.func.name}`.toLowerCase();
    const i = keywords.findIndex((k) => hay.includes(k));
    return i === -1 ? keywords.length : i;
  };
  return [...items].sort((a, b) => {
    const ap = a.packageName === primary ? 0 : 1;
    const bp = b.packageName === primary ? 0 : 1;
    if (ap !== bp) return ap - bp;
    const ar = priorityRank(a);
    const br = priorityRank(b);
    if (ar !== br) return ar - br;
    const pkg = a.packageName.localeCompare(b.packageName);
    if (pkg !== 0) return pkg;
    return a.name.localeCompare(b.name);
  });
}

/** MIME-ish key used to carry a node type name through HTML5 drag/drop.
 *  The canvas drop handler reads this and adds the corresponding node at
 *  the drop point. */
export const FF_DRAG_MIME = 'application/x-funcflow-node';

export interface FunctionBrowserCallbacks {
  onFunctionDoubleClick: (funcInfo: FuncInfo) => void;
  onBuiltinNodeDoubleClick: (nodeTypeName: string) => void;
  /** Double-click (or Enter) on a file in the Files tree → add an OpenFile node. */
  onFileDoubleClick: (file: DG.FileInfo) => void;
  /** Files chosen via the Files-tab Upload button → add Uploaded File nodes. */
  onLocalFilesPicked: (files: File[]) => void;
}

/** Connection display name for a query func, used to group queries in the
 *  Queries pane (`friendlyName` falls back to `name`). */
export function queryConnectionName(f: FuncInfo): string {
  if (!(f.func instanceof DG.DataQuery)) return '';
  try {
    const conn = (f.func as DG.DataQuery).connection as DG.DataConnection | null;
    return (conn?.friendlyName ?? conn?.name ?? 'Other') || 'Other';
  } catch {
    return 'Other';
  }
}

/** Persisted UI state for the browser (the group-by mode). Section
 *  expand/collapse state is NOT kept here — the platform accordion
 *  (`ui.accordion(key)`) persists it itself in `localStorage['Accordion:<key>']`,
 *  keyed by pane name. */
interface BrowserState {
  groupBy: GroupByMode;
}
const LS_KEY = 'funcflow.browser.v1';
/** Persistence key of the toolbox accordion; the nested per-connection
 *  accordion inside the Queries pane uses `<key>.queries`. */
export const TOOLBOX_ACCORDION_KEY = 'funcflow.toolbox';
/** Persistence key of the top tab strip (Files / Queries / Workflows /
 *  Favorites) — DG.TabControl remembers the selected tab itself. */
export const TOOLBOX_TABS_KEY = 'funcflow.toolbox.tab';
/** The top-strip tab names, in order. */
export const TOOLBOX_TABS = ['Files', 'Queries', 'Workflows', 'Favorites'] as const;
const VALID_MODES: GroupByMode[] = ['category', 'role', 'tags', 'package'];

/** Case- AND whitespace-insensitive substring match — so "openfile" matches
 *  "Open File" and "table output" matches "TableOutput". */
export function nameMatchesQuery(text: string, query: string): boolean {
  if (!query) return true;
  const t = (text || '').toLowerCase();
  const q = query.toLowerCase();
  return t.includes(q) || t.replace(/\s+/g, '').includes(q.replace(/\s+/g, ''));
}

/** Whether a function matches a search query. Matches the display name AND the
 *  raw function name / friendlyName — the toolbox shows the friendly name (e.g.
 *  "Open File") but users search by the name they know (e.g. "OpenFile"), and
 *  vice versa — all whitespace-insensitive; plus description, tags, role, pkg. */
export function funcMatchesSearch(f: FuncInfo, query: string): boolean {
  if (!query) return true;
  const q = query.toLowerCase();
  return nameMatchesQuery(f.name, q) ||
    nameMatchesQuery(f.func.name || '', q) ||
    nameMatchesQuery(f.func.friendlyName || '', q) ||
    (f.func.description || '').toLowerCase().includes(q) ||
    f.tags.some((t) => t.toLowerCase().includes(q)) ||
    (f.role || '').toLowerCase().includes(q) ||
    f.packageName.toLowerCase().includes(q);
}

/** Left sidebar: searchable, groupable function catalog */
/** Display labels for the group-by modes (the "by: …" catalog-header button). */
const GROUP_BY_LABELS: Record<GroupByMode, string> = {
  category: 'what it does',
  role: 'role',
  tags: 'tags',
  package: 'package',
};

/** FA5 icon per toolbox section header (muted grey — recognition, not decor).
 *  Covers the "what it does" categories and the fixed panes; other group-by
 *  modes produce keys outside this map and simply render without an icon. */
const SECTION_ICONS: Record<string, string> = {
  'Cheminformatics': 'flask',
  'Bioinformatics': 'dna',
  'Data Sources': 'database',
  'Combine Tables': 'layer-group',
  'Transform Tables': 'random',
  'Column Operations': 'columns',
  'Compute Values': 'calculator',
  'Visualize': 'chart-bar',
  'Other': 'ellipsis-h',
  'Viewers': 'chart-pie',
  'Widgets': 'puzzle-piece',
  'Inputs': 'sign-in-alt',
  'Outputs': 'sign-out-alt',
  'Constants': 'hashtag',
  'Utilities': 'wrench',
  'Debug': 'bug',
};

/** A muted section-header icon (fixed 16px box so labels align). `name` is a
 *  bare FA5 name (light weight) or a full `fas fa-…` override. */
function sectionIcon(name: string): HTMLElement {
  const i = document.createElement('i');
  const cls = name.includes(' ') ? name : `fal fa-${name}`;
  i.className = `grok-icon ${cls} funcflow-section-icon`;
  return i;
}

export class FunctionBrowser {
  root: HTMLElement;
  private searchInput!: HTMLInputElement;
  private groupByBtn!: HTMLElement;
  private treeContainer!: HTMLElement;
  private groupBy: GroupByMode = 'category';
  /** The toolbox accordion, rebuilt on every render. A keyed `ui.accordion`
   *  restores pane states from localStorage at construction and persists them
   *  on header clicks — no manual bookkeeping. Exposed for tests. */
  accordion: DG.Accordion | null = null;
  /** The nested per-connection accordion inside the Queries pane (built lazily
   *  when that pane first expands). Exposed for tests. */
  queriesAccordion: DG.Accordion | null = null;
  private callbacks: FunctionBrowserCallbacks;
  /** The Files tree (KNIME-style file browser). Built lazily once and reused
   *  across renders so its expanded state and scroll position survive a search. */
  private filesTreeRoot?: HTMLElement;
  /** The top tab strip (Files / Queries / Workflows / Favorites). Exposed for tests. */
  topTabs: DG.TabControl | null = null;
  private queriesTabContent!: HTMLElement;
  private workflowsTabContent!: HTMLElement;
  private favoritesTabContent!: HTMLElement;
  private favoritesUnsub: (() => void) | null = null;

  constructor(callbacks: FunctionBrowserCallbacks) {
    this.callbacks = callbacks;
    this.loadState();
    this.root = this.buildUI();
    // A star toggled anywhere updates every visible star and the Favorites tab
    // in place — no full re-render, so tree scroll positions survive.
    this.favoritesUnsub = onFavoritesChanged(() => {
      this.syncStars();
      this.renderFavoritesTab();
    });
  }

  /** Unhooks the global favorites listener. Call when the hosting view closes. */
  destroy(): void {
    this.favoritesUnsub?.();
    this.favoritesUnsub = null;
  }

  // ---------- persistence (localStorage) ----------

  private loadState(): void {
    try {
      const raw = localStorage.getItem(LS_KEY);
      if (!raw) return;
      const s = JSON.parse(raw) as Partial<BrowserState>;
      if (s.groupBy && VALID_MODES.includes(s.groupBy)) this.groupBy = s.groupBy;
    } catch {/* corrupt/blocked storage — fall back to defaults */}
  }

  private saveState(): void {
    try {
      const state: BrowserState = {groupBy: this.groupBy};
      localStorage.setItem(LS_KEY, JSON.stringify(state));
    } catch {/* storage blocked/full — non-fatal */}
  }

  private buildUI(): HTMLElement {
    // Search bar — filters everything below it: the collection tabs
    // (queries / workflows / favorites, with match counts) AND the catalog.
    this.searchInput = document.createElement('input');
    this.searchInput.type = 'text';
    this.searchInput.placeholder = 'Search functions, queries, flows…';
    this.searchInput.className = 'funcflow-search-input';
    setTid(this.searchInput, 'browser-search');
    this.searchInput.addEventListener('input', () => this.render());

    // Clear (✕) affordance at the right edge of the search box — shown only
    // while there's text to clear.
    const clearBtn = document.createElement('span');
    clearBtn.className = 'funcflow-search-clear';
    clearBtn.innerHTML = '&times;';
    ui.tooltip.bind(clearBtn, 'Clear search');
    setTid(clearBtn, 'browser-search-clear');
    const syncClear = (): void => {clearBtn.style.display = this.searchInput.value ? 'flex' : 'none';};
    clearBtn.addEventListener('click', () => {
      this.searchInput.value = '';
      this.searchInput.dispatchEvent(new Event('input', {bubbles: true}));
      this.searchInput.focus();
      syncClear();
    });
    this.searchInput.addEventListener('input', syncClear);
    syncClear();
    const searchWrap = ui.div([this.searchInput, clearBtn], 'funcflow-search-wrap');

    // Group-by: a compact "by: <mode> ▾" text-button in the catalog zone
    // header — it scopes ONLY the categories accordion, so it sits on that
    // zone rather than next to the global search. Opens a popup menu.
    this.groupByBtn = setTid(ui.div([], 'funcflow-groupby-btn'), 'browser-groupby');
    this.syncGroupByLabel();
    ui.tooltip.bind(this.groupByBtn, 'How the function list below is organized');
    this.groupByBtn.addEventListener('click', () => {
      DG.Menu.popup()
        .items(Object.values(GROUP_BY_LABELS), (label) => {
          const mode = (Object.keys(GROUP_BY_LABELS) as GroupByMode[])
            .find((k) => GROUP_BY_LABELS[k] === label)!;
          this.setGroupBy(mode);
        }, {radioGroup: 'ff-groupby', isChecked: (label) => GROUP_BY_LABELS[this.groupBy] === label})
        .show();
    });
    const zoneLabel = ui.div([], 'funcflow-zone-label');
    zoneLabel.textContent = 'Functions';
    const catalogHeader = setTid(
      ui.div([zoneLabel, this.groupByBtn], 'funcflow-catalog-header'), 'browser-catalog-header');

    this.treeContainer = setTid(ui.div([], 'funcflow-tree-container'), 'browser-tree');

    // Zones top-to-bottom: global search, the collection tabs, the labeled
    // function catalog, and (appended by the view) the Suggestions strip.
    const container = ui.divV([
      searchWrap,
      this.buildTopTabs(),
      catalogHeader,
      this.treeContainer,
    ], 'funcflow-browser');

    return setTid(container, 'browser');
  }

  /** Switches the accordion grouping mode (the popup menu + tests use this). */
  setGroupBy(mode: GroupByMode): void {
    this.groupBy = mode;
    this.saveState();
    this.syncGroupByLabel();
    this.render();
  }

  private syncGroupByLabel(): void {
    this.groupByBtn.textContent = `by: ${GROUP_BY_LABELS[this.groupBy]} ▾`;
  }

  /** The top strip: a platform tab control with the item *collections* — Files,
   *  Queries, Workflows (saved flows), and Favorites (starred nodes) — leaving
   *  the accordion below to the function categories. The selected tab persists
   *  via the TabControl key; the Files tree builds lazily on first show. */
  private buildTopTabs(): HTMLElement {
    const tabs = DG.TabControl.create(false, TOOLBOX_TABS_KEY);
    this.topTabs = tabs;
    this.queriesTabContent = setTid(ui.div([], 'funcflow-tab-content'), 'browser-queries');
    this.workflowsTabContent = setTid(ui.div([], 'funcflow-tab-content'), 'browser-workflows');
    this.favoritesTabContent = setTid(ui.div([], 'funcflow-tab-content'), 'browser-favorites');
    const filesContent = setTid(ui.div([], 'funcflow-tab-content'), 'browser-files');

    const panes: Array<{name: string; content: () => HTMLElement; tip: string}> = [
      {name: 'Files', tip: 'Browse data files. Double-click or drag a file onto the canvas to load it.',
        content: () => {
          if (filesContent.childElementCount === 0) {
            const text = document.createElement('span');
            text.className = 'funcflow-tab-hint-text';
            // Ends by pointing at the open-local-file button beside it; short
            // enough to share the row with it in two lines.
            text.textContent = 'Double-click or drag a file to load it — or open a local one:';
            const hint = ui.div([text, this.buildUploadButton()], 'funcflow-tab-hint');
            filesContent.appendChild(hint);
          }
          filesContent.appendChild(this.getFilesTreeRoot());
          return filesContent;
        }},
      {name: 'Queries', tip: 'Database queries, grouped by data connection. Double-click or drag to add.',
        content: () => this.queriesTabContent},
      {name: 'Workflows', tip: 'Saved flows — reuse them as nodes in this flow. Double-click or drag to add.',
        content: () => this.workflowsTabContent},
      {name: 'Favorites', tip: 'Nodes you starred. Hover any node in the toolbox and click its ★ to pin it here.',
        content: () => this.favoritesTabContent},
    ];
    for (const p of panes) {
      const pane = tabs.addPane(p.name, p.content);
      setTid(pane.header, 'browser-tab', p.name);
      pane.header.dataset.tab = p.name;
      // Wrap the header text so a squeezed tab truncates its LABEL with an
      // ellipsis instead of clipping the search badge appended after it
      // (narrow toolbox / zoomed-in screens).
      const label = document.createElement('span');
      label.className = 'funcflow-tab-label';
      while (pane.header.firstChild) label.appendChild(pane.header.firstChild);
      pane.header.appendChild(label);
      ui.tooltip.bind(pane.header, p.tip);
    }

    const host = ui.div([tabs.root], 'funcflow-top-tabs');
    return setTid(host, 'browser-tabs');
  }

  /** Programmatically activates one of the top tabs (used by the guide). */
  showTab(name: (typeof TOOLBOX_TABS)[number]): void {
    if (this.topTabs) this.topTabs.currentPane = this.topTabs.getPane(name);
  }

  /** Domain sections floated to the very top of the function list (right after
   *  the Queries pane) — the science a chemist/biologist reaches for first. */
  private static readonly DOMAIN_CATEGORIES = ['Cheminformatics', 'Bioinformatics'];

  render(): void {
    // While a query is active the Suggestions strip (canvas-contextual, not
    // search results) mutes via CSS so its items aren't mistaken for matches.
    this.root.dataset.searching = String(!!this.searchInput.value);
    this.treeContainer.innerHTML = '';
    // The platform accordion: a keyed one restores/persists pane states itself.
    // While a search is active every matching section is forced open — build an
    // UNKEYED accordion then, so the forced states don't overwrite the user's.
    const acc = ui.accordion(this.searchInput.value ? null : TOOLBOX_ACCORDION_KEY);
    this.accordion = acc;
    this.queriesAccordion = null;

    // DG function nodes — queries, workflows, widget-producers, and (already-
    // excluded) viewer functions don't appear here; queries and workflows live
    // in the top tabs, widgets/viewers in their dedicated panes.
    const funcs = this.filterBySearch(getRegisteredFuncs()
      .filter((f) => !(f.func instanceof DG.DataQuery) && !funcOutputsWidget(f) && !isWorkflowFunc(f.func)));
    const grouped = this.groupFunctions(funcs);
    const renderCategory = (category: string): void => {
      const raw = grouped[category];
      if (!raw || raw.length === 0) return;
      const items = FunctionBrowser.DOMAIN_CATEGORIES.includes(category) ?
        orderDomainSection(raw, category) : raw;
      this.addSection(acc, category, () => {
        const content = ui.div([], 'funcflow-section-content');
        for (const info of items) content.appendChild(this.createFuncItem(info));
        return content;
      }, {count: items.length});
    };

    // Section order: the Cheminformatics / Bioinformatics domain sections (the
    // science a chemist/biologist reaches for first), the task categories
    // (Data Sources → … → Compute Values), Viewers, Widgets, the
    // building-block built-ins (Inputs / Outputs / Constants / Utilities),
    // then the Other catch-all, and Debug last. Files, Queries, and Workflows
    // live in the top tabs, not here.
    for (const domain of FunctionBrowser.DOMAIN_CATEGORIES) renderCategory(domain);

    const sortedKeys = this.orderGroupKeys(Object.keys(grouped))
      .filter((k) => !FunctionBrowser.DOMAIN_CATEGORIES.includes(k));
    for (const category of sortedKeys.filter((k) => k !== 'Other')) renderCategory(category);

    this.renderViewersSection(acc);
    this.renderWidgetsSection(acc);
    this.renderBuiltinNodes(acc, ['Inputs', 'Outputs', 'Constants', 'Utilities']);
    if (sortedKeys.includes('Other')) renderCategory('Other');
    this.renderBuiltinNodes(acc, ['Debug']);

    acc.end();
    this.treeContainer.appendChild(acc.root);
    this.renderTabs();
  }

  /** Refreshes the search-dependent top-tab contents (Files is a persistent
   *  tree and refreshes itself). */
  private renderTabs(): void {
    this.renderQueriesTab();
    this.renderWorkflowsTab();
    this.renderFavoritesTab();
    this.updateTabBadges();
  }

  /** During a search each tab header shows how many of its rows match (blue
   *  count badge; a 0-match tab dims) — matches hiding behind an inactive tab
   *  stay discoverable without auto-switching tabs. Files isn't part of the
   *  search and stays neutral. Cleared when the query clears. */
  private updateTabBadges(): void {
    if (!this.topTabs) return;
    const query = this.searchInput.value.toLowerCase();
    const counts: Record<string, number> = {
      'Queries': getRegisteredFuncs()
        .filter((f) => f.func instanceof DG.DataQuery && funcMatchesSearch(f, query)).length,
      'Workflows': getRegisteredFuncs()
        .filter((f) => isWorkflowFunc(f.func) && funcMatchesSearch(f, query)).length,
      'Favorites': getFavorites().filter((e) => nameMatchesQuery(e.label, query)).length,
    };
    for (const name of Object.keys(counts)) {
      let header: HTMLElement | null = null;
      try {
        header = this.topTabs.getPane(name)?.header ?? null;
      } catch {/* pane missing — nothing to badge */}
      if (!header) continue;
      let badge = header.querySelector<HTMLElement>('.funcflow-tab-badge');
      // A 0-match tab just dims — a "0" badge adds noise and crowds the strip.
      if (!query || counts[name] === 0) {
        badge?.remove();
        header.classList.toggle('funcflow-tab-dim', !!query && counts[name] === 0);
        continue;
      }
      if (!badge) {
        badge = document.createElement('span');
        badge.className = 'funcflow-tab-badge';
        header.appendChild(badge);
      }
      badge.textContent = String(counts[name]);
      header.classList.remove('funcflow-tab-dim');
    }
  }

  /** Adds one toolbox section as an accordion pane. Content builds lazily on
   *  first expand; on a keyed accordion the remembered state wins over the
   *  `expanded` default. While a search is active the pane is forced open (the
   *  accordion is unkeyed then, so nothing persists). Pane headers keep the
   *  `data-section` attribute and `ff-browser-section-<title>` test id the
   *  guide and tests address sections by. */
  private addSection(acc: DG.Accordion, title: string, getContent: () => HTMLElement, opts: {
    expanded?: boolean; count?: number; tooltip?: string; tidParts?: Array<string | number>;
  } = {}): DG.AccordionPane {
    const pane = opts.count !== undefined ?
      acc.addCountPane(title, getContent, () => opts.count!, opts.expanded ?? false, null, false) :
      acc.addPane(title, getContent, opts.expanded ?? false, null, false);
    if (this.searchInput.value) pane.expanded = true;
    const header = pane.root.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
    if (header) {
      header.dataset.section = title;
      setTid(header, 'browser-section', title);
      if (SECTION_ICONS[title]) header.insertBefore(sectionIcon(SECTION_ICONS[title]), header.firstChild);
      if (opts.tooltip) ui.tooltip.bind(header, opts.tooltip);
    }
    if (opts.tidParts) setTid(pane.root, ...opts.tidParts);
    return pane;
  }

  /** In "what it does" mode, order groups by the curated task sequence
   *  (Data Sources first); in every other mode fall back to alphabetical. */
  private orderGroupKeys(keys: string[]): string[] {
    if (this.groupBy !== 'category') return keys.sort();
    const idx = (k: string): number => {
      const i = (FUNC_CATEGORIES as readonly string[]).indexOf(k);
      return i === -1 ? FUNC_CATEGORIES.length : i;
    };
    return keys.sort((a, b) => idx(a) - idx(b) || a.localeCompare(b));
  }

  /** The KNIME-style Files browser — built once, reused across renders so its
   *  expanded folders and scroll position survive a search keystroke. */
  private getFilesTreeRoot(): HTMLElement {
    if (!this.filesTreeRoot) {
      const tree = getFilesBrowser(
        () => {},                                                  // selection: no-op
        (f) => this.callbacks.onFileDoubleClick(f),                // dbl-click / Enter → OpenFile node
        'funcflow.files.v1');
      this.filesTreeRoot = tree.root;
      // The tab pane owns the height now — drop the tree's standalone cap.
      this.filesTreeRoot.style.maxHeight = '';
    }
    return this.filesTreeRoot;
  }

  /** The Files-tab "Upload" chip: opens the OS file picker (filtered to the
   *  table formats the platform can read) and hands the picked files to the
   *  view, which registers the bytes and adds an Uploaded File node per file —
   *  the same pipeline as dropping files onto the canvas. Lives inside the
   *  hint row, so it costs no extra vertical space. */
  private buildUploadButton(): HTMLElement {
    const input = document.createElement('input');
    input.type = 'file';
    input.multiple = true;
    input.style.display = 'none';
    input.addEventListener('change', () => {
      const files = Array.from(input.files ?? []);
      if (files.length > 0) this.callbacks.onLocalFilesPicked(files);
      input.value = ''; // re-picking the same file must fire change again
    });
    // input.click() bubbles back up to the chip — don't re-trigger the picker.
    input.addEventListener('click', (ev) => ev.stopPropagation());
    // The platform's "Open local file" affordance is a folder icon — reuse
    // that vocabulary here, icon-only (the tooltip carries the explanation).
    const icon = document.createElement('i');
    icon.className = 'grok-icon fal fa-folder-open';
    const btn = ui.div([icon, input], 'funcflow-upload-btn');
    ui.tooltip.bind(btn, 'Open a local file — it becomes a node on the canvas ' +
      'and is stored with the flow when you save.');
    btn.addEventListener('click', () => {
      // Resolved at click time so handlers from lazily-loaded packages count.
      input.accept = supportedUploadExtensions().map((e) => `.${e}`).join(',');
      input.click();
    });
    return setTid(btn, 'browser-upload');
  }

  /** A muted hint shown when a top tab has nothing to list. */
  private static emptyNote(text: string): HTMLElement {
    const el = ui.div([], 'funcflow-tab-empty');
    el.textContent = text;
    return el;
  }

  private renderQueriesTab(): void {
    // Group every registered query by its connection (friendlyName ?? name).
    const queries = getRegisteredFuncs().filter((f) => f.func instanceof DG.DataQuery);
    const query = this.searchInput.value.toLowerCase();
    const matching = query ? queries.filter((f) => funcMatchesSearch(f, query)) : queries;
    this.queriesTabContent.innerHTML = '';
    this.queriesAccordion = null;
    if (matching.length === 0) {
      this.queriesTabContent.appendChild(FunctionBrowser.emptyNote(
        query ? 'No queries match the search.' : 'No database queries available.'));
      return;
    }

    const byConn = new Map<string, FuncInfo[]>();
    for (const f of matching) {
      const conn = queryConnectionName(f);
      let arr = byConn.get(conn);
      if (!arr) byConn.set(conn, arr = []);
      arr.push(f);
    }
    this.queriesTabContent.appendChild(this.buildQueriesContent(byConn));
  }

  private renderWorkflowsTab(): void {
    const flows = getRegisteredFuncs().filter((f) => isWorkflowFunc(f.func));
    const query = this.searchInput.value.toLowerCase();
    const matching = query ? flows.filter((f) => funcMatchesSearch(f, query)) : flows;
    matching.sort((a, b) => a.name.localeCompare(b.name));
    this.workflowsTabContent.innerHTML = '';
    if (matching.length === 0) {
      this.workflowsTabContent.appendChild(FunctionBrowser.emptyNote(
        query ? 'No workflows match the search.' :
          'No saved flows yet. Save a flow and it appears here, ready to reuse as a node.'));
      return;
    }
    const content = ui.div([], 'funcflow-section-content');
    for (const info of matching) content.appendChild(this.createFuncItem(info));
    this.workflowsTabContent.appendChild(content);
  }

  private renderFavoritesTab(): void {
    const query = this.searchInput.value.toLowerCase();
    const favs = getFavorites().filter((e) => nameMatchesQuery(e.label, query));
    this.favoritesTabContent.innerHTML = '';
    if (favs.length === 0) {
      this.favoritesTabContent.appendChild(FunctionBrowser.emptyNote(
        query ? 'No favorites match the search.' :
          'Nothing starred yet. Hover any node in the toolbox and click its ★ to pin it here.'));
      return;
    }
    const content = ui.div([], 'funcflow-section-content');
    for (const e of favs) {
      // A favorite may be a DG function (rich FuncInfo available) or a
      // builtin/viewer type — either way the type name alone creates the node.
      const info = getRegisteredFuncs().find((f) => f.nodeTypeName === e.type);
      const item = this.makeToolboxItem(e.label, e.type);
      item.dataset.testid = tid('browser-fav-item', e.type);
      if (info) item.dataset.func = info.func.name;
      ui.tooltip.bind(item, info?.func.description || `${e.label}. Double-click or drag to add.`);
      item.addEventListener('dblclick', () => {
        if (info) this.callbacks.onFunctionDoubleClick(info);
        else this.callbacks.onBuiltinNodeDoubleClick(e.type);
      });
      this.makeItemDraggable(item, e.type);
      content.appendChild(item);
    }
    this.favoritesTabContent.appendChild(content);
  }

  /** The Queries pane content: a nested accordion with one pane per connection
   *  (its own persistence key — pane names are the connection names). */
  private buildQueriesContent(byConn: Map<string, FuncInfo[]>): HTMLElement {
    const hasSearch = !!this.searchInput.value;
    const inner = ui.accordion(hasSearch ? null : `${TOOLBOX_ACCORDION_KEY}.queries`);
    this.queriesAccordion = inner;
    for (const conn of [...byConn.keys()].sort((a, b) => a.localeCompare(b))) {
      const items = byConn.get(conn)!;
      items.sort((a, b) => a.name.localeCompare(b.name));
      const pane = inner.addCountPane(conn, () => {
        const subContent = ui.div([], 'funcflow-section-content');
        for (const info of items)
          subContent.appendChild(this.createFuncItem(info));
        return subContent;
      }, () => items.length, false, null, false);
      if (hasSearch) pane.expanded = true;
      pane.root.dataset.queryConn = conn;
      setTid(pane.root, 'browser-query-conn', conn);
      const header = pane.root.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
      if (header) {
        // A uniform glyph down the list says "these are all connections" —
        // together with the tray background it keeps the Queries tab from
        // reading as more function categories.
        header.insertBefore(sectionIcon('database'), header.firstChild);
        ui.tooltip.bind(header, `Queries from the “${conn}” connection`);
      }
    }
    inner.end();
    return ui.div([inner.root], 'funcflow-section-content funcflow-subsections');
  }

  /** The Viewers pane — manually-built viewer nodes (core charts first, then
   *  discovered package viewers). Each adds a `Viewers/<label>` node. */
  private renderViewersSection(acc: DG.Accordion): void {
    const query = this.searchInput.value.toLowerCase();
    // "chart"/"plot"/"graph" are what a scientist actually types when hunting
    // for a visualization — any of them surfaces the whole Viewers pane.
    const synonymHit = query.length >= 3 &&
      ['chart', 'charts', 'plot', 'plots', 'graph', 'graphs', 'viewer', 'viewers']
        .some((w) => w.startsWith(query) || query.startsWith(w));
    const types = query && !synonymHit ?
      VIEWER_NODE_TYPES.filter((v) => nameMatchesQuery(v.label, query)) : VIEWER_NODE_TYPES;
    if (types.length === 0) return;

    this.addSection(acc, 'Viewers', () => {
      const content = ui.div([], 'funcflow-section-content');
      for (const vt of types) {
        const item = this.makeToolboxItem(vt.label, vt.nodeTypeName);
        item.dataset.testid = tid('browser-item', vt.nodeTypeName);
        ui.tooltip.bind(item, `Add a ${vt.label} viewer. Wire a table in, then run. Double-click or drag.`);
        item.addEventListener('dblclick', () => this.callbacks.onBuiltinNodeDoubleClick(vt.nodeTypeName));
        this.makeItemDraggable(item, vt.nodeTypeName);
        content.appendChild(item);
      }
      return content;
    }, {
      count: types.length,
      tooltip: 'Charts and viewers. Wire a table into one and run to see it in the preview panel.',
      tidParts: ['browser-viewers'],
    });
  }

  /** The Widgets pane — functions that produce a widget (info panels, search
   *  widgets, …), kept out of the categories. */
  private renderWidgetsSection(acc: DG.Accordion): void {
    const widgets = getRegisteredFuncs().filter(funcOutputsWidget);
    const query = this.searchInput.value.toLowerCase();
    const matching = query ? widgets.filter((f) => funcMatchesSearch(f, query)) : widgets;
    if (matching.length === 0) return;
    matching.sort((a, b) => a.name.localeCompare(b.name));

    this.addSection(acc, 'Widgets', () => {
      const content = ui.div([], 'funcflow-section-content');
      for (const info of matching)
        content.appendChild(this.createFuncItem(info));
      return content;
    }, {
      count: matching.length,
      tooltip: 'Functions that produce a widget. Double-click or drag to add.',
      tidParts: ['browser-widgets'],
    });
  }

  /** Renders the built-in node sections named in `titles` (lets the caller
   *  interleave them with other panes — Debug goes last in the toolbox). */
  private renderBuiltinNodes(acc: DG.Accordion, titles: string[]): void {
    const query = this.searchInput.value.toLowerCase();

    // Input nodes
    const inputNodes = [
      {name: 'Table Input', type: 'Inputs/Table Input', desc: 'Dataframe input parameter'},
      {name: 'Column Input', type: 'Inputs/Column Input', desc: 'Single column from a table'},
      {name: 'Column List Input', type: 'Inputs/Column List Input', desc: 'Multiple columns from a table'},
      {name: 'String Input', type: 'Inputs/String Input', desc: 'Text input parameter'},
      {name: 'Number Input', type: 'Inputs/Number Input', desc: 'Floating-point number input'},
      {name: 'Int Input', type: 'Inputs/Int Input', desc: 'Integer number input'},
      {name: 'Boolean Input', type: 'Inputs/Boolean Input', desc: 'True/false toggle input'},
      {name: 'DateTime Input', type: 'Inputs/DateTime Input', desc: 'Date and time input'},
      {name: 'File Input', type: 'Inputs/File Input', desc: 'File upload input'},
      {name: 'Map Input', type: 'Inputs/Map Input', desc: 'Key-value map input'},
      {name: 'Dynamic Input', type: 'Inputs/Dynamic Input', desc: 'Dynamically-typed input'},
      {name: 'String List Input', type: 'Inputs/String List Input', desc: 'List of strings input'},
      {name: 'Blob Input', type: 'Inputs/Blob Input', desc: 'Binary data input'},
    ];
    const outputNodes = [
      {name: 'Table Output', type: 'Outputs/Table Output', desc: 'Marks a dataframe as script output'},
      {name: 'Value Output', type: 'Outputs/Value Output', desc: 'Marks a value as script output (configurable type)'},
    ];
    const utilityNodes = [
      {name: 'Select Column', type: 'Utilities/Select Column', desc: 'Gets a column from a table by name'},
      {name: 'Select Columns', type: 'Utilities/Select Columns', desc: 'Gets multiple columns by names (comma-separated)'},
      {name: 'Select Table', type: 'Utilities/Select Table', desc: 'Gets an open table by name via grok.shell.tableByName()'},
      {name: 'Add Table View', type: 'Utilities/Add Table View', desc: 'Opens a table in a new view via grok.shell.addTableView()'},
      {name: 'Log', type: 'Utilities/Log', desc: 'Logs a value to the browser console'},
      {name: 'Info', type: 'Utilities/Info', desc: 'Shows an info balloon via grok.shell.info()'},
      {name: 'Warning', type: 'Utilities/Warning', desc: 'Shows a warning balloon via grok.shell.warning()'},
      {name: 'ToString', type: 'Utilities/ToString', desc: 'Converts any value to a string'},
      {name: 'FromJSON', type: 'Utilities/FromJSON', desc: 'Parses a JSON string into an object'},
      {name: 'ToJSON', type: 'Utilities/ToJSON', desc: 'Serializes a value to a JSON string'},
    ];
    const constantNodes = [
      {name: 'String', type: 'Constants/String', desc: 'A constant text value'},
      {name: 'Int', type: 'Constants/Int', desc: 'A constant integer value'},
      {name: 'Double', type: 'Constants/Double', desc: 'A constant floating-point value'},
      {name: 'Boolean', type: 'Constants/Boolean', desc: 'A constant true/false value'},
      {name: 'List', type: 'Constants/List', desc: 'A constant list of comma-separated values'},
    ];
    const debugNodes = [
      {name: 'Breakpoint', type: 'Debug/Breakpoint', desc: 'Pauses execution in debug mode until Continue is clicked'},
    ];

    // All built-in sections start collapsed: they are building-blocks a scientist
    // reaches for deliberately, not the first thing to scan. The data-function
    // categories below (Data Sources first) are what's expanded-ready instead.
    const sections: {title: string; nodes: {name: string; type: string; desc?: string}[]; tip?: string}[] = [
      {title: 'Inputs', nodes: inputNodes, tip: 'Script input parameters (become //input: lines)'},
      {title: 'Outputs', nodes: outputNodes, tip: 'Script output parameters (become //output: lines)'},
      {title: 'Constants', nodes: constantNodes, tip: 'Constant literal values'},
      {title: 'Utilities', nodes: utilityNodes, tip: 'Helper operations (logging, type conversion, etc.)'},
      {title: 'Debug', nodes: debugNodes, tip: 'Debugging and execution control nodes'},
    ];

    for (const section of sections.filter((s) => titles.includes(s.title))) {
      const filtered = query ?
        section.nodes.filter((n) => nameMatchesQuery(n.name, query)) :
        section.nodes;
      if (filtered.length === 0) continue;

      this.addSection(acc, section.title, () => {
        const content = ui.div([], 'funcflow-section-content');
        for (const node of filtered) {
          const item = this.makeToolboxItem(node.name, node.type); // type e.g. "Inputs/Table Input"
          item.dataset.testid = tid('browser-item', node.type);
          const tip = node.desc ?
            `${node.desc}. Double-click or drag to add` :
            `Double-click or drag to add ${node.name}`;
          ui.tooltip.bind(item, tip);
          item.addEventListener('dblclick', () => {
            this.callbacks.onBuiltinNodeDoubleClick(node.type);
          });
          this.makeItemDraggable(item, node.type);
          content.appendChild(item);
        }
        return content;
      }, {tooltip: section.tip, count: filtered.length});
    }
  }

  /** Sets a star icon's visual state (outline vs. filled gold). */
  private static applyStarState(star: HTMLElement, fav: boolean): void {
    star.className = `funcflow-item-star grok-icon ${fav ? 'fas' : 'fal'} fa-star` +
      (fav ? ' funcflow-item-star-active' : '');
  }

  /** The base toolbox row every catalog item shares: an ellipsized label plus a
   *  trailing ★ that stars the node type into the Favorites tab. Sets
   *  `data-node-type-name` (used by node creation, drag, and star sync). */
  private makeToolboxItem(label: string, typeName: string): HTMLElement {
    const labelEl = ui.div([], 'funcflow-item-label');
    labelEl.textContent = label;
    const star = document.createElement('i');
    FunctionBrowser.applyStarState(star, isFavorite(typeName));
    ui.tooltip.bind(star, () => isFavorite(typeName) ? 'Remove from Favorites' : 'Add to Favorites');
    star.addEventListener('click', (ev) => {
      ev.stopPropagation();
      toggleFavorite({type: typeName, label});
    });
    star.addEventListener('dblclick', (ev) => ev.stopPropagation());
    setTid(star, 'browser-item-star', typeName);
    const item = ui.div([labelEl, star], 'funcflow-func-item');
    item.dataset.nodeTypeName = typeName;
    return item;
  }

  /** Repaints every visible star after a favorites change. */
  private syncStars(): void {
    for (const item of Array.from(this.root.querySelectorAll<HTMLElement>('.funcflow-func-item[data-node-type-name]'))) {
      const star = item.querySelector<HTMLElement>('.funcflow-item-star');
      if (star) FunctionBrowser.applyStarState(star, isFavorite(item.dataset.nodeTypeName!));
    }
  }

  /** One draggable, double-clickable catalog row for a DG function. Shared by
   *  the category sections and the Queries pane. */
  private createFuncItem(info: FuncInfo): HTMLElement {
    const item = this.makeToolboxItem(info.name, info.nodeTypeName);
    item.dataset.testid = tid('browser-item', info.nodeTypeName);
    item.dataset.func = info.func.name;                // the DG function it adds
    if (info.packageName) item.dataset.package = info.packageName;
    let tip = info.func.description || info.name;
    if (info.packageName)
      tip += ` (${info.packageName})`;
    ui.tooltip.bind(item, tip);
    item.addEventListener('dblclick', () => {
      this.callbacks.onFunctionDoubleClick(info);
    });
    this.makeItemDraggable(item, info.nodeTypeName);
    return item;
  }

  /** Wire HTML5 drag on a browser item so dropping it on the canvas creates
   *  the matching node at the drop point. */
  private makeItemDraggable(item: HTMLElement, typeName: string): void {
    item.draggable = true;
    item.addEventListener('dragstart', (ev) => {
      if (!ev.dataTransfer) return;
      ev.dataTransfer.setData(FF_DRAG_MIME, typeName);
      ev.dataTransfer.setData('text/plain', typeName);
      ev.dataTransfer.effectAllowed = 'copy';
      item.classList.add('funcflow-func-item-dragging');
    });
    item.addEventListener('dragend', () => {
      item.classList.remove('funcflow-func-item-dragging');
    });
  }

  private filterBySearch(funcs: FuncInfo[]): FuncInfo[] {
    const query = this.searchInput.value.toLowerCase();
    if (!query) return funcs;
    return funcs.filter((f) => funcMatchesSearch(f, query));
  }

  private groupFunctions(funcs: FuncInfo[]): Record<string, FuncInfo[]> {
    const groups: Record<string, FuncInfo[]> = {};

    for (const f of funcs) {
      let keys: string[];
      // Saved flows get their own 'Workflows' section in EVERY grouping — their
      // role/tags/package say nothing useful (a script from whatever package
      // the user happened to save it in).
      if (isWorkflowFunc(f.func))
        keys = ['Workflows'];
      else {
        switch (this.groupBy) {
        case 'role':
          keys = [f.role || 'Uncategorized'];
          break;
        case 'tags':
          keys = f.tags.length > 0 ? f.tags : ['Untagged'];
          break;
        case 'package':
          keys = [f.packageName || 'Core'];
          break;
        case 'category':
          keys = [categorizeFunc(f.func, f.role, f.packageName)];
          break;
        default:
          keys = ['Other'];
        }
      }
      for (const key of keys) {
        if (!groups[key]) groups[key] = [];
        groups[key].push(f);
      }
    }

    return groups;
  }
}
