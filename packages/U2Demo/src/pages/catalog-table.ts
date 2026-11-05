/* The catalog table over 50 made-up design campaigns: every choice it offers (sort by a double click,
   groups with counts, the columns menu with a submenu of fields, row actions by access, a context
   menu), its first load, a refresh that brings a campaign made elsewhere, and its empty and error
   states. The choices are kept in this browser (`persist`). */
import {
  signal, computed, AsyncValue,
  divV, divH, span, h1, dot, button, notify,
  CatalogTable, ButtonGroup, SortOrder,
} from '@datagrok-libraries/u2';
import type {CatalogColumn, CatalogGrouping, CatalogAccess, Action, BadgeVariant, ReadonlySignal}
  from '@datagrok-libraries/u2';
import {readout} from './common';

interface Campaign {
  code: string;
  name: string;
  template: string;
  author: string;
  modifiedBy: string;
  created: string;
  modified: number;
  molecules: number;
  status: string;
  fields: {[name: string]: unknown};
}

const TEMPLATES: {[key: string]: {name: string, fields: string[]}} = {
  KHL: {name: 'Kinase hit-to-lead', fields: ['Target', 'Chemist', 'Potency goal (nM)', 'Kickoff']},
  GPCR: {name: 'GPCR series', fields: ['Receptor', 'Chemist', 'Selectivity window (fold)']},
  FRAG: {name: 'Fragment expansion', fields: ['Fragment library', 'Hit rate (%)', 'Covalent']},
  SKP: {name: 'Quick sketchpad', fields: []},
};
const PEOPLE = ['Admin', 'skalkin', 'maria.dolotova', 'greg.landrum', 'pavel', 'alex.aprm'];
const STATUSES = ['In Progress', 'In Progress', 'In Progress', 'Submitted', 'Completed', 'On hold', 'Archived'];
const TARGETS = ['EGFR', 'JAK2', 'CDK9', 'IRAK4', 'PI3Kδ', 'BTK', 'ALK', 'LRRK2'];
const TOPICS = ['selectivity push', 'back-up series', 'quick wins', 'macrocycles', 'degrader warheads',
  'brain-penetrant', 'covalent follow-up', 'resistance mutants', 'linking', 'fragment merges'];
const LONG = ' — kinome-wide selectivity, hERG margin and solubility';

function campaigns(count: number, seed = 7): Campaign[] {
  let s = seed;
  const next = () => (s = (s * 16807) % 2147483647) / 2147483647;
  const pick = <X>(xs: X[]): X => xs[Math.floor(next() * xs.length)];
  const numbers: {[key: string]: number} = {};
  return Array.from({length: count}, (_, i) => {
    const template = pick(Object.keys(TEMPLATES));
    const n = numbers[template] = (numbers[template] ?? 0) + 1;
    const fields: {[name: string]: unknown} = {};
    for (const f of TEMPLATES[template].fields) {
      if (next() < 0.2)
        continue;
      fields[f] = f === 'Kickoff' ? new Date(2026, Math.floor(next() * 10), 1 + Math.floor(next() * 27)) :
        f === 'Covalent' ? next() < 0.5 : f.includes('(') ? Math.round(next() * 100) : f === 'Chemist' ? pick(PEOPLE) :
          pick(TARGETS);
    }
    return {
      code: `${template}-${n}`,
      name: `${pick(TARGETS)} ${pick(TOPICS)}${i % 9 === 4 ? LONG : ''}`,
      template,
      author: pick(PEOPLE),
      modifiedBy: pick(PEOPLE),
      created: `2026/${1 + Math.floor(next() * 10)}/${1 + Math.floor(next() * 28)}`,
      modified: Math.floor(next() * 1e6),
      molecules: Math.floor(next() * next() * 1500),
      status: i % 13 === 5 ? '' : pick(STATUSES),
      fields,
    };
  });
}

const STATUS_DOT: {[status: string]: BadgeVariant} = {'In Progress': 'accent', 'Submitted': 'success',
  'Completed': 'success', 'On hold': 'warning'};

const fieldText = (v: unknown): string => v instanceof Date ? v.toLocaleDateString() :
  typeof v === 'boolean' ? (v ? 'Yes' : 'No') : v === undefined || v === null ? '' : String(v);

function columnsOf(items: readonly Campaign[]): CatalogColumn<Campaign>[] {
  const fields = [...new Set(items.flatMap((c) => Object.keys(c.fields)))].sort();
  const owners = (f: string) => Object.entries(TEMPLATES).filter(([, t]) => t.fields.includes(f))
    .map(([, t]) => t.name).join(', ');
  return [
    {name: 'Name', value: (c) => c.name, pinned: true, link: 'Continue Campaign', maxWidth: '320px'},
    {name: 'Code', value: (c) => c.code},
    {name: 'Created', value: (c) => c.created},
    {name: 'Author', value: (c) => c.author},
    {name: 'Last Modified by', value: (c) => c.modifiedBy},
    {name: 'Last Modified', value: (c) => new Date(1761000000000 + c.modified * 1000), type: 'date'},
    {name: 'Molecules', value: (c) => c.molecules, align: 'right'},
    {name: 'Status', value: (c) => c.status,
      render: (c) => c.status === '' ? '' : divH([dot(STATUS_DOT[c.status] ?? 'default'), span(c.status)], {gap: 's'})},
    ...fields.map((f): CatalogColumn<Campaign> => ({name: `campaignFields.${f}`, header: f,
      value: (c) => c.fields[f] instanceof Date ? c.fields[f] : fieldText(c.fields[f]),
      render: (c) => fieldText(c.fields[f]), menuGroup: 'Campaign Fields', menuDetail: owners(f)})),
  ];
}

function groupingsOf(items: readonly Campaign[]): CatalogGrouping<Campaign>[] {
  const fields = [...new Set(items.flatMap((c) => Object.keys(c.fields)))].sort();
  return [
    {name: 'Template', key: (c) => c.template, groupLabel: (key) => TEMPLATES[key]?.name ?? key},
    {name: 'Status', key: (c) => c.status},
    {name: 'Author', key: (c) => c.author},
    {name: 'Last Modified User', key: (c) => c.modifiedBy},
    ...fields.map((f): CatalogGrouping<Campaign> => ({name: `campaignFields.${f}`, label: f,
      key: (c) => fieldText(c.fields[f]), menuGroup: 'Campaign Fields'})),
  ];
}

export function catalogTablePage(): HTMLElement {
  const mode = signal<'data' | 'empty' | 'error'>('data');
  let all = campaigns(50);
  let added = 0;
  const wait = (ms: number) => new Promise((resolve) => setTimeout(resolve, ms));
  // what the columns and groupings are made of: the last answer, kept while a newer one comes
  const known = signal<readonly Campaign[]>([]);
  const source = new AsyncValue<readonly Campaign[]>(async () => {
    await wait(600);
    if (mode.peek() === 'error')
      throw new Error('the campaigns folder cannot be read');
    return known.value = mode.peek() === 'empty' ? [] : all;
  }, {debounceMs: 0});

  // the edit rights arrive later, all at once, as one lookup per group would bring them
  const rights = signal<'pending' | 'known'>('pending');
  void wait(900).then(() => rights.value = 'known');
  const access = new Map<string, ReadonlySignal<CatalogAccess>>();
  const accessOf = (c: Campaign): ReadonlySignal<CatalogAccess> => {
    let a = access.get(c.code);
    if (a === undefined) {
      a = computed(() => rights.value === 'pending' ? 'pending' : c.author === 'alex.aprm' ? 'view' : 'edit');
      access.set(c.code, a);
    }
    return a;
  };

  const table: CatalogTable<Campaign> = new CatalogTable<Campaign>({
    name: 'campaigns',
    label: 'Campaigns',
    items: source,
    keyOf: (c) => c.code,
    columns: computed(() => columnsOf(known.value)),
    groupings: computed(() => groupingsOf(known.value)),
    defaultColumns: ['Code', 'Created', 'Author', 'Last Modified by', 'Molecules', 'Status'],
    defaultOrder: (a, b) => b.modified - a.modified || a.code.localeCompare(b.code),
    emptyGroupLabel: 'No Information',
    actionsVisible: 'always',
    access: accessOf,
    viewOnly: {icon: 'eye', tooltip: 'View only: you can open this campaign but not change it'},
    actions: (c): Action[] => [
      {text: 'Manage campaign permissions', icon: 'users', requires: 'edit',
        run: () => void notify.info(`Permissions of ${c.code} would open here`)},
      {text: 'Delete campaign', icon: 'trash-alt', requires: 'edit', destructive: true, run: () => {
        all = all.filter((x) => x.code !== c.code);
        void source.refresh().then(() => notify.info(`Campaign ${c.code} deleted`));
      }},
    ],
    onActivate: (c) => notify.info(`${c.code} would open here`),
    empty: 'No campaigns yet. Start the first one under New Campaign.',
    footer: (shown) => {
      const viewOnly = shown.filter((c) => accessOf(c).value === 'view').length;
      return `${shown.length} campaigns${viewOnly > 0 ? ` · ${viewOnly} view only` : ''}`;
    },
    persist: 'u2demo-campaigns',
  });

  const chrome = table.chrome({columns: 'Edit visible columns in campaigns table',
    group: (g) => `Group Campaigns. Current: ${g}`, refresh: 'Refresh campaigns table'});
  const show = (m: 'data' | 'empty' | 'error') => {
    mode.value = m;
    void source.refresh();
  };
  const states = new ButtonGroup({toggle: 'single', items: [
    {id: 'data', text: 'Campaigns', onClick: () => show('data')},
    {id: 'empty', text: 'None', onClick: () => show('empty')},
    {id: 'error', text: 'Unreadable', onClick: () => show('error')},
  ]});
  states.selected.value = ['data'];
  const elsewhere = () => {
    added++;
    all = [{...campaigns(1, 1000 + added)[0], code: `SKP-${100 + added}`, name: `Made elsewhere ${added}`,
      modified: 2e6 + added}, ...all];
    void table.refresh();
  };
  const sort = computed(() => {
    const s = SortOrder.parse(table.state.value.sort);
    return s === null ? '(last modified)' : `${s.column} ${s.descending ? 'descending' : 'ascending'}`;
  });

  return divV([
    span('A catalog of records: double-click a header to sort (again for descending, a third time for the ' +
      'default order), group from the layer icon, choose columns from the eye; right-click a row for its ' +
      'actions. The choices are kept in this browser.', 'u2demo-hint'),
    divH([h1('Continue Campaign'), ...chrome], 'u2demo-catalog-heading'),
    table,
    divH([span('Show'), states, span('a campaign made elsewhere appears on refresh:'),
      button('Add one elsewhere', elsewhere)], 'u2demo-row'),
    readout('sort', sort),
    readout('groupBy', computed(() => table.state.value.groupBy ?? 'None')),
    readout('columns', computed(() => table.state.value.columns.join(', '))),
  ], 'u2demo-page u2demo-page-wide');
}
