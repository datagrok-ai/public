/* eslint-disable max-len */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {
  FuncInfo, getRegisteredFuncs, isWorkflowFunc, ensureFuncNodeType, loadQueryFuncs, VIEWER_NODE_TYPES,
} from '../rete/node-factory';
import {BUILTIN_SECTIONS} from '../rete/builtin-catalog';
import {tid, setTid} from '../utils/test-ids';
import {getFilesBrowser} from '../utils/files-browser-tree';
import {SpacePicker} from '../ui/space-picker';
import {categorizeBySignature, domainCategory} from '../types/type-map';
import {getFavorites, isFavorite, toggleFavorite, onFavoritesChanged} from './favorites';
import {supportedUploadExtensions} from '../utils/uploaded-files';

export function funcOutputsWidget(f: FuncInfo): boolean {
  try {
    return f.func.outputs.some((o) => String(o.propertyType) === 'widget');
  } catch {
    return false;
  }
}

export type GroupByMode = 'category' | 'role' | 'tags' | 'package';

/** Task categories in pipeline order — used both to bucket and to order the tree. */
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

/** Bucket a function by what it does — domain package wins over the signature category. */
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

/** Flagship package whose functions lead each domain section. */
export const DOMAIN_PRIMARY_PACKAGE: Record<string, string> = {
  Cheminformatics: 'Chem',
  Bioinformatics: 'Bio',
};

/** Lowercased substring matches against the function name; lower index sorts higher, unlisted last. */
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

/** Drag-payload key carrying a node type name through HTML5 drag/drop. */
export const FF_DRAG_MIME = 'application/x-funcflow-node';

export interface FunctionBrowserCallbacks {
  onFunctionDoubleClick: (funcInfo: FuncInfo) => void;
  onBuiltinNodeDoubleClick: (nodeTypeName: string) => void;
  onFileDoubleClick: (file: DG.FileInfo) => void;
  onLocalFilesPicked: (files: File[]) => void;
}

/** Queries-pane grouping key; empty when the query carries no connection (such queries are skipped). */
export function queryConnectionName(f: FuncInfo): string {
  if (!(f.func instanceof DG.DataQuery)) return '';
  try {
    const conn = (f.func as DG.DataQuery).connection as DG.DataConnection | null;
    return conn ? ((conn.friendlyName ?? conn.name) || '') : '';
  } catch {
    return '';
  }
}

/** Only the group-by mode — pane expand state is persisted by the keyed accordion itself. */
interface BrowserState {
  groupBy: GroupByMode;
}
const LS_KEY = 'funcflow.browser.v1';
/** Persistence key of the toolbox accordion; the Queries pane nests `<key>.queries`. */
export const TOOLBOX_ACCORDION_KEY = 'funcflow.toolbox';
/** Persistence key of the top tab strip — DG.TabControl remembers the selected tab itself. */
export const TOOLBOX_TABS_KEY = 'funcflow.toolbox.tab';
export const TOOLBOX_TABS = ['Files', 'Spaces', 'Queries', 'Workflows', 'Favorites'] as const;

/** Same glyphs the platform Browse tree uses for these concepts (browse_panel_tree.dart). */
export const TOOLBOX_TAB_ICONS: Record<(typeof TOOLBOX_TABS)[number], string> = {
  'Files': 'folder',
  'Spaces': 'brackets-curly',
  'Queries': 'database',
  'Workflows': 'sitemap',
  'Favorites': 'star',
};
const VALID_MODES: GroupByMode[] = ['category', 'role', 'tags', 'package'];

/** Case- and whitespace-insensitive substring match. */
export function nameMatchesQuery(text: string, query: string): boolean {
  if (!query) return true;
  const t = (text || '').toLowerCase();
  const q = query.toLowerCase();
  return t.includes(q) || t.replace(/\s+/g, '').includes(q.replace(/\s+/g, ''));
}

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

const GROUP_BY_LABELS: Record<GroupByMode, string> = {
  category: 'what it does',
  role: 'role',
  tags: 'tags',
  package: 'package',
};

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

/** `name` is a bare FA5 name (light weight) or a full `fas fa-…` override. */
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
  /** Rebuilt on every render; a keyed accordion persists pane states itself. */
  accordion: DG.Accordion | null = null;
  queriesAccordion: DG.Accordion | null = null;
  private callbacks: FunctionBrowserCallbacks;
  private filesTreeRoot?: HTMLElement;
  topTabs: DG.TabControl | null = null;
  private queriesTabContent!: HTMLElement;
  private workflowsTabContent!: HTMLElement;
  private favoritesTabContent!: HTMLElement;
  private spacesTabContent!: HTMLElement;
  spacePicker: SpacePicker | null = null;
  private favoritesUnsub: (() => void) | null = null;
  private queryFuncs: FuncInfo[] | null = null;
  private queryFuncsRequested = false;

  constructor(callbacks: FunctionBrowserCallbacks) {
    this.callbacks = callbacks;
    this.loadState();
    this.root = this.buildUI();
    // Stars update in place — no full re-render, so tree scroll positions survive.
    this.favoritesUnsub = onFavoritesChanged(() => {
      this.syncStars();
      this.renderFavoritesTab();
    });
  }

  destroy(): void {
    this.favoritesUnsub?.();
    this.favoritesUnsub = null;
  }

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
    this.searchInput = document.createElement('input');
    this.searchInput.type = 'text';
    this.searchInput.placeholder = 'Search functions, queries, flows…';
    this.searchInput.className = 'funcflow-search-input';
    setTid(this.searchInput, 'browser-search');
    this.searchInput.addEventListener('input', () => this.render());

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

    const container = ui.divV([
      searchWrap,
      this.buildTopTabs(),
      catalogHeader,
      this.treeContainer,
    ], 'funcflow-browser');

    return setTid(container, 'browser');
  }

  setGroupBy(mode: GroupByMode): void {
    this.groupBy = mode;
    this.saveState();
    this.syncGroupByLabel();
    this.render();
  }

  private syncGroupByLabel(): void {
    this.groupByBtn.textContent = `by: ${GROUP_BY_LABELS[this.groupBy]} ▾`;
  }

  private buildTopTabs(): HTMLElement {
    const tabs = DG.TabControl.create(false, TOOLBOX_TABS_KEY);
    this.topTabs = tabs;
    this.queriesTabContent = setTid(ui.div([], 'funcflow-tab-content'), 'browser-queries');
    this.workflowsTabContent = setTid(ui.div([], 'funcflow-tab-content'), 'browser-workflows');
    this.favoritesTabContent = setTid(ui.div([], 'funcflow-tab-content'), 'browser-favorites');
    this.spacesTabContent = setTid(ui.div([], 'funcflow-tab-content'), 'browser-spaces');
    const filesContent = setTid(ui.div([], 'funcflow-tab-content'), 'browser-files');

    const panes: Array<{name: (typeof TOOLBOX_TABS)[number]; content: () => HTMLElement; tip: string}> = [
      {name: 'Files', tip: 'Files — browse data files. Double-click or drag a file onto the canvas to load it.',
        content: () => {
          if (filesContent.childElementCount === 0) {
            const text = document.createElement('span');
            text.className = 'funcflow-tab-hint-text';
            text.textContent = 'Double-click or drag a file to load it — or open a local one:';
            const hint = ui.div([text, this.buildUploadButton()], 'funcflow-tab-hint');
            filesContent.appendChild(hint);
          }
          filesContent.appendChild(this.getFilesTreeRoot());
          return filesContent;
        }},
      {name: 'Spaces', tip: 'Spaces — shared hierarchical folders. Double-click a flow, query, or file inside ' +
        'a space to use it here.',
        content: () => this.getSpacesTabContent()},
      {name: 'Queries', tip: 'Queries — database queries, grouped by data connection. Double-click or drag to add.',
        content: () => this.queriesTabContent},
      {name: 'Workflows', tip: 'Workflows — saved flows, reusable as nodes in this flow. Double-click or drag to add.',
        content: () => this.workflowsTabContent},
      {name: 'Favorites', tip: 'Favorites — nodes you starred. Hover any node in the toolbox and click its ★ ' +
        'to pin it here.',
        content: () => this.favoritesTabContent},
    ];
    for (const p of panes) {
      const pane = tabs.addPane(p.name, p.content);
      setTid(pane.header, 'browser-tab', p.name);
      pane.header.dataset.tab = p.name;
      // The label span is kept so the search badge appended after it never clips.
      const label = document.createElement('span');
      label.className = 'funcflow-tab-label';
      while (pane.header.firstChild) pane.header.firstChild.remove();
      const icon = document.createElement('i');
      icon.className = `grok-icon fal fa-${TOOLBOX_TAB_ICONS[p.name]} funcflow-tab-icon`;
      label.appendChild(icon);
      pane.header.appendChild(label);
      pane.header.setAttribute('aria-label', p.name);
      ui.tooltip.bind(pane.header, p.tip);
    }

    const host = ui.div([tabs.root], 'funcflow-top-tabs');
    return setTid(host, 'browser-tabs');
  }

  /** Built once on first activation and reused so expanded spaces survive a search. */
  private getSpacesTabContent(): HTMLElement {
    if (this.spacesTabContent.childElementCount === 0) {
      const hint = ui.div([], 'funcflow-tab-hint');
      hint.textContent = 'Browse spaces — double-click or drag a flow, query, or file to add it.';
      this.spacesTabContent.appendChild(hint);
      const host = ui.div([ui.loader()]);
      this.spacesTabContent.appendChild(host);
      void SpacePicker.create({
        showContent: true,
        onEntityActivated: (e) => this.activateSpaceEntity(e),
      }).then((picker) => {
        this.spacePicker = picker;
        host.replaceChildren(picker.root);
      });
    }
    return this.spacesTabContent;
  }

  private activateSpaceEntity(e: DG.Entity): void {
    if (e instanceof DG.FileInfo) {
      if (e.isFile) this.callbacks.onFileDoubleClick(e);
      return;
    }
    if (!(e instanceof DG.Func)) return;
    const info = getRegisteredFuncs().find((f) => f.func.id === e.id);
    if (info)
      this.callbacks.onFunctionDoubleClick(info);
    else
      this.callbacks.onBuiltinNodeDoubleClick(ensureFuncNodeType(e));
  }

  showTab(name: (typeof TOOLBOX_TABS)[number]): void {
    if (this.topTabs) this.topTabs.currentPane = this.topTabs.getPane(name);
  }

  private static readonly DOMAIN_CATEGORIES = ['Cheminformatics', 'Bioinformatics'];

  render(): void {
    this.root.dataset.searching = String(!!this.searchInput.value);
    this.treeContainer.innerHTML = '';
    // Search forces matching panes open — build an UNKEYED accordion then, so forced states never overwrite the user's.
    const acc = ui.accordion(this.searchInput.value ? null : TOOLBOX_ACCORDION_KEY);
    this.accordion = acc;
    this.queriesAccordion = null;

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

  private renderTabs(): void {
    this.renderQueriesTab();
    this.renderWorkflowsTab();
    this.renderFavoritesTab();
    this.updateTabBadges();
  }

  private updateTabBadges(): void {
    if (!this.topTabs) return;
    const query = this.searchInput.value.toLowerCase();
    const counts: Record<string, number> = {
      'Queries': (this.queryFuncs ?? []).filter((f) => funcMatchesSearch(f, query)).length,
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

  /** One toolbox section pane; content builds lazily on first expand, and an active search forces it open. */
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

  private orderGroupKeys(keys: string[]): string[] {
    if (this.groupBy !== 'category') return keys.sort();
    const idx = (k: string): number => {
      const i = (FUNC_CATEGORIES as readonly string[]).indexOf(k);
      return i === -1 ? FUNC_CATEGORIES.length : i;
    };
    return keys.sort((a, b) => idx(a) - idx(b) || a.localeCompare(b));
  }

  /** Built once and reused so expanded folders and scroll survive a search keystroke. */
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

  private static emptyNote(text: string): HTMLElement {
    const el = ui.div([], 'funcflow-tab-empty');
    el.textContent = text;
    return el;
  }

  private renderQueriesTab(): void {
    // The query list comes from the server (loadQueryFuncs) — the registry's DG.Func.find scan misses queries.
    if (this.queryFuncs === null) {
      this.queriesTabContent.innerHTML = '';
      this.queriesAccordion = null;
      this.queriesTabContent.appendChild(ui.loader());
      if (!this.queryFuncsRequested) {
        this.queryFuncsRequested = true;
        loadQueryFuncs().catch((e) => {
          console.warn('FuncFlow: query catalog load failed', e);
          return [] as FuncInfo[];
        }).then((funcs) => {
          this.queryFuncs = funcs;
          this.renderQueriesTab();
          this.updateTabBadges();
        });
      }
      return;
    }
    const query = this.searchInput.value.toLowerCase();
    const matching = query ? this.queryFuncs.filter((f) => funcMatchesSearch(f, query)) : this.queryFuncs;
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
      if (!conn) continue; // no connection object — skipped by design
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
        header.insertBefore(sectionIcon('database'), header.firstChild);
        ui.tooltip.bind(header, `Queries from the “${conn}” connection`);
      }
    }
    inner.end();
    return ui.div([inner.root], 'funcflow-section-content funcflow-subsections');
  }

  private renderViewersSection(acc: DG.Accordion): void {
    const query = this.searchInput.value.toLowerCase();
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

  private renderBuiltinNodes(acc: DG.Accordion, titles: string[]): void {
    const query = this.searchInput.value.toLowerCase();

    for (const section of BUILTIN_SECTIONS.filter((s) => titles.includes(s.title))) {
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

  private static applyStarState(star: HTMLElement, fav: boolean): void {
    star.className = `funcflow-item-star grok-icon ${fav ? 'fas' : 'fal'} fa-star` +
      (fav ? ' funcflow-item-star-active' : '');
  }

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

  private syncStars(): void {
    for (const item of Array.from(this.root.querySelectorAll<HTMLElement>('.funcflow-func-item[data-node-type-name]'))) {
      const star = item.querySelector<HTMLElement>('.funcflow-item-star');
      if (star) FunctionBrowser.applyStarState(star, isFavorite(item.dataset.nodeTypeName!));
    }
  }

  private createFuncItem(info: FuncInfo): HTMLElement {
    const item = this.makeToolboxItem(info.name, info.nodeTypeName);
    item.dataset.testid = tid('browser-item', info.nodeTypeName);
    item.dataset.func = info.func.name;
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
      // Saved flows get their own 'Workflows' section in EVERY grouping — their role/tags/package say nothing useful.
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
