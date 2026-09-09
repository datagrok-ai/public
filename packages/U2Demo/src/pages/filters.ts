/* The filter builder and the query input over the three schema sources the platform layer
   offers: a DataFrame (client-side, evaluated into `df.filter`), an entity type (the smart-filter
   string a dapi `filter()` takes) and a domain table (the condition tree a domain query takes).
   The active tab and its query ride the app URL as `?tab=<id>&q=<query>`, so a filter is a link. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {
  signal, computed, Scope,
  divV, divH, span, h3, button,
  TabStrip, TextInput, AsyncView, FilterBuilder, FilterQueryInput, Filters,
} from '@datagrok-libraries/u2';
import type {Signal, FilterGroup, FilterSchema} from '@datagrok-libraries/u2';
import {FilterSchemas, toBitSet, chip, EntityChip, viewers} from '@datagrok-libraries/u2/src/dg/index.js';
import {demoOrders} from '../package';
import {readout} from './common';

type FilterTab = 'dataframe' | 'entity' | 'domain';
const TABS: FilterTab[] = ['dataframe', 'entity', 'domain'];

/** The page's URL state: the active tab and one query per tab. The inbound path is written here
 * before the page builds and read back as each pane's first tree — applied as soon as the pane's
 * schema is ready — so the filter survives a visit to another sub-demo, and a copied link opens
 * on the same tab with the same rows. */
const activeTab = signal<FilterTab>('dataframe');
const queries: Record<FilterTab, Signal<string>> = {dataframe: signal(''), entity: signal(''), domain: signal('')};

/** `base?tab=<id>&q=<query>` — the tab omitted on the first one, the query when empty; encoded
 * the way `Filters.queryPath` encodes. */
export function filterPath(base: string): string {
  const tab = activeTab.value;
  const parts: string[] = [];
  if (tab !== 'dataframe')
    parts.push(`tab=${tab}`);
  const query = queries[tab].value;
  if (query !== '')
    parts.push(`q=${encodeURIComponent(query)}`);
  return parts.length === 0 ? base : `${base}?${parts.join('&')}`;
}

/** Resets the state to the path's `tab` and `q` (a path or a `location.search`); false where the
 * path carries neither. */
export function readFilterPath(path: string): boolean {
  const at = path.indexOf('?');
  const params = new URLSearchParams(at < 0 ? '' : path.slice(at + 1).split('#')[0]);
  const tab = params.get('tab');
  activeTab.value = TABS.find((t) => t === tab) ?? 'dataframe';
  for (const t of TABS)
    queries[t].value = '';
  queries[activeTab.peek()].value = params.get('q') ?? '';
  return params.has('tab') || params.has('q');
}

/** The two editors over one tree, with the query readout — the shape every tab repeats. */
function editors(schema: FilterSchema, tree: Signal<FilterGroup>, name: string, placeholder: string,
  target: 'dataframe' | 'domain'): HTMLElement[] {
  const fb = new FilterBuilder({label: 'Criteria', schema, bind: tree, name, target, showQuery: true});
  const q = new FilterQueryInput({label: 'Query', schema, bind: tree, name: `${name}Query`, target, placeholder});
  return [q.root, fb.root];
}

/** A pane's tree: the tab's query where it parses, else empty; the tab's query follows the tree
 * from then on. */
function restored(scope: Scope, tab: FilterTab, schema: FilterSchema, target: 'dataframe' | 'domain'):
  {tree: Signal<FilterGroup>, problem: string} {
  const initial = Filters.parse(queries[tab].peek(), schema, target);
  const tree = signal(initial.problems.length === 0 ? initial.root : Filters.group('and'));
  scope.effect(() => queries[tab].value = Filters.format(tree.value));
  return {tree, problem: initial.problems[0]?.message ?? ''};
}

/** A pane that needs the server: the schema loads once, the content builds inside the view's own
 * scope (a lazy tab or a click has no ambient one), and where the platform cannot answer (local
 * mode, no such table) the failure is the pane's content, with Retry. */
function asyncPane(load: () => Promise<FilterSchema>, render: (schema: FilterSchema) => HTMLElement[]):
  AsyncView<FilterSchema> {
  const view = AsyncView.owned(async () => [await load()],
    ([schema]) => divV(render(schema), 'u2demo-filters-pane'));
  view.root.classList.add('u2demo-filters-pane');
  view.refresh();
  return view;
}

function dataFramePane(scope: Scope): HTMLElement {
  const df = demoOrders(90);
  const schema = FilterSchemas.forDataFrame(df);
  const {tree, problem} = restored(scope, 'dataframe', schema, 'dataframe');

  const grid = viewers.grid(df);
  grid.root.classList.add('u2demo-filters-grid');
  const rows = signal(`${df.rowCount} of ${df.rowCount}`);
  const message = signal(problem);
  const apply = async (): Promise<void> => {
    const root = tree.peek();
    const problems = Filters.validate(root, schema, 'dataframe');
    message.value = problems[0]?.message ?? '';
    if (problems.length > 0)
      return;
    if (Filters.count(root) === 0)
      df.filter.setAll(true);
    else
      df.filter.copyFrom(await toBitSet(df, root));
    rows.value = `${df.filter.trueCount} of ${df.rowCount}`;
  };
  if (Filters.count(tree.peek()) > 0)
    void apply();

  return divV([
    span('Six demo orders, filtered right here in the browser: describe the rows you want, and ' +
      'Apply narrows the table below.', 'u2demo-hint'),
    ...editors(schema, tree, 'orders', 'city = "Basel" and total > 500', 'dataframe'),
    divH([button('Apply', () => void apply(), {primary: true}),
      button('Clear', () => tree.value = Filters.group('and'))], 'u2demo-row'),
    readout('rows', rows),
    readout('problem', computed(() => message.value || '(none)')),
    grid.root,
  ], 'u2demo-filters-pane');
}

function entityTypePane(): HTMLElement {
  return asyncPane(() => FilterSchemas.forEntityType('User'), (schema) => {
    const scope = Scope.ambient!;
    const {tree} = restored(scope, 'entity', schema, 'domain');
    const count = signal('(not run)');
    const picked = divH([], 'u2demo-row');
    const chips: EntityChip[] = [];
    scope.own(() => {
      for (const c of chips)
        c.dispose();
    });
    const run = async (): Promise<void> => {
      const users = await grok.dapi.users.filter(Filters.format(tree.peek())).list({pageSize: 20});
      for (const c of chips)
        c.dispose();
      chips.length = 0;
      picked.textContent = '';
      for (const user of users) {
        const c = chip(user);
        chips.push(c);
        picked.append(c.root);
      }
      count.value = String(users.length);
    };
    const start = () => void run().catch((e) => count.value = String(e));
    if (Filters.count(tree.peek()) > 0)
      start();
    return [
      span('Find the platform\'s users by what is known about them — their login, when they joined, ' +
        'the groups they belong to.', 'u2demo-hint'),
      ...editors(schema, tree, 'users', 'login like "adm" and joined > -1y', 'domain'),
      divH([button('Run', start, {primary: true})], 'u2demo-row'),
      readout('users', count),
      picked,
    ];
  }).root;
}

function domainTablePane(): HTMLElement {
  const scope = Scope.ambient!;
  const address = new TextInput({label: 'Table', placeholder: 'schema.table', name: 'domainAddress'});
  const host = divV([], 'u2demo-filters-pane');
  let pane: AsyncView<FilterSchema> | undefined;
  scope.own(() => pane?.dispose());
  const open = (): void => {
    const a = address.value.peek();
    pane?.dispose();
    pane = asyncPane(() => FilterSchemas.forDomainTable(a), (schema) => {
      const {tree} = restored(Scope.ambient!, 'domain', schema, 'domain');
      const count = signal('(not run)');
      const run = async (): Promise<void> => {
        const root = tree.peek();
        const q = grok.dapi.domains.table(a).query();
        // u2's platform-free tree type and the API's are the same shape, spelled apart
        const where = Filters.toDomainTree(root) as DG.DomainFilter;
        count.value = String(await (Filters.count(root) === 0 ? q : q.where(where)).count());
      };
      const start = () => void run().catch((e) => count.value = String(e));
      if (Filters.count(tree.peek()) > 0)
        start();
      return [
        ...editors(schema, tree, 'domain', 'name like "a" and created_on > -1m', 'domain'),
        divH([button('Count', start, {primary: true})], 'u2demo-row'),
        readout('rows', count),
      ];
    });
    host.replaceChildren(pane.root);
  };
  grok.meta.coreLocationOf('User').then((at) => {
    if (at == null || address.value.peek() !== '')
      return;
    address.value.value = `${at.schema}.${at.table}`;
    // a link into this tab names its rows: open the table it was made over
    if (queries.domain.peek() !== '')
      open();
  }, () => {});
  return divV([
    span('Filter the rows of any table the platform knows: name the table, open it, and Count ' +
      'tells how many rows match.', 'u2demo-hint'),
    divH([address.root, button('Open', open)], 'u2demo-row'),
    host,
  ], 'u2demo-filters-pane');
}

export function filtersPage(): HTMLElement {
  const scope = Scope.ambient!;
  const tabs = new TabStrip({tabs: [
    {id: 'dataframe', label: 'DataFrame', content: dataFramePane(scope)},
    // a lazy pane builds on the first click, outside any ambient scope: run it under the page's
    {id: 'entity', label: 'Entity type', content: () => Scope.runWith(scope, entityTypePane)},
    {id: 'domain', label: 'Domain table', content: () => Scope.runWith(scope, domainTablePane)},
  ]});
  tabs.name = 'filterTabs';
  tabs.activeTab.value = activeTab.peek();
  scope.effect(() => activeTab.value = TABS.find((t) => t === tabs.activeTab.value) ?? 'dataframe');
  return divV([
    span('One immutable FilterGroup, two editors: the builder (rows of property, operator, ' +
      'value under one and/or) and the query input (the grammar with completion). Each tab ' +
      'binds both to a schema from a different source.', 'u2demo-hint'),
    h3('Filter builder'),
    tabs,
  ], 'u2demo-page');
}
