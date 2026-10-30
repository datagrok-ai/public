/* The campaigns table of an app's landing: the u2 CatalogTable over the campaign.json files the user may view, with
   the columns, groupings and row actions of the HitTriage apps. Hit Design and PeptiHit keep their choices under the
   HD* localStorage keys, Hit Triage and PepTriage under the HT* ones, in the formats the landings always used. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {AsyncValue, CatalogTable, SortOrder, confirmDialog, computed, divH, dot, plural, signal, span}
  from '@datagrok-libraries/u2';
import type {Action, BadgeVariant, CatalogAccess, CatalogColumn, CatalogGrouping, CatalogPersistence, CatalogState,
  ReadonlySignal} from '@datagrok-libraries/u2';
import {_package} from '../../package';
import * as api from '../../package-api';
import {AppName, HitDesignCampaign, HitTriageCampaign} from '../types';
import {CampaignGrouping, CampaignJsonName, DefaultCampaignTableInfoGetters, HDCampaignsCollapsedLSKey,
  HDCampaignsGroupingLSKey, HDCampaignsTableSortingLSKey, HDCampaignTableColumnsLSKey, HDDefaultCampaignTableColumns,
  HTCampaignsCollapsedLSKey, HTCampaignsGroupingLSKey, HTCampaignsTableSortingLSKey, HTCampaignTableColumnsLSKey,
  HTDefaultCampaignTableColumns, HTDefaultCampaignTableInfoGetters, i18n, LastModifiedByColumn} from '../consts';
import {checkEditPermissions, MembershipCache} from '../utils';
import {defaultPermissions, PermissionsEditor} from './permissions-editor';

export type AnyCampaign = HitDesignCampaign | HitTriageCampaign;

/** Hit Design and PeptiHit ('HD'), Hit Triage and PepTriage ('HT'): each pair shares its table settings. */
export type CampaignFamily = 'HD' | 'HT';

interface Family {
  keys: {grouping: string, columns: string, sorting: string, collapsed: string};
  getters: {[column: string]: (c: AnyCampaign) => string};
  defaultColumns: string[];
  counts: string[];
  /** The last save of each campaign is in hitdesign.update_logs. */
  saveLog: boolean;
}

const FAMILIES: {[f in CampaignFamily]: Family} = {
  HD: {keys: {grouping: HDCampaignsGroupingLSKey, columns: HDCampaignTableColumnsLSKey,
    sorting: HDCampaignsTableSortingLSKey, collapsed: HDCampaignsCollapsedLSKey},
  getters: DefaultCampaignTableInfoGetters as Family['getters'], defaultColumns: HDDefaultCampaignTableColumns,
  counts: ['Molecules'], saveLog: true},
  HT: {keys: {grouping: HTCampaignsGroupingLSKey, columns: HTCampaignTableColumnsLSKey,
    sorting: HTCampaignsTableSortingLSKey, collapsed: HTCampaignsCollapsedLSKey},
  getters: HTDefaultCampaignTableInfoGetters as unknown as Family['getters'],
  defaultColumns: HTDefaultCampaignTableColumns, counts: ['Total', 'Selected'], saveLog: false},
};

const LEGACY_GROUP_KEY = 'Accordion:Hit-Design-campaigns-group-';
const LAST_MODIFIED = 'Last Modified';
const FIELD = 'campaignFields.';
const CODES = new Intl.Collator(undefined, {numeric: true});

export interface CampaignCatalogOptions {
  app: AppName;
  family: CampaignFamily;
  /** Campaigns deleted in this page: a share listing may still show their folders for a while. */
  deleted: string[];
  /** The name link: opens the campaign. */
  open: (code: string) => unknown;
  /** Deletes the campaign's files and drops it from the page's list. */
  remove: (code: string) => Promise<void>;
  /** How the balloons name a campaign: Hit Design by its code, Hit Triage by its name. */
  label: (c: AnyCampaign) => string;
  /** Every list read, the statuses of a Submit dialog's suggestions among what it feeds. */
  onLoaded?: (campaigns: readonly AnyCampaign[]) => void;
  /** The table's host is narrow (under 820px): a name is cut at 200px rather than 320px. */
  narrow?: ReadonlySignal<boolean>;
}

export interface CampaignCatalog {
  readonly table: CatalogTable<AnyCampaign>;
  /** The heading's icons: columns, grouping and refresh, with the landing's tooltips. */
  chrome(): HTMLElement[];
  /** Reads the list again: from the page's copy, or with `fresh` from the server. The rows stay meanwhile. */
  reload(fresh?: boolean): Promise<void>;
}

/** The table and its data: a list read through the page's campaign cache (Refresh reads the server again), the
 * order of the last save, and the user's edit rights, looked up once per group for the whole list. Built where it
 * is owned: its sources go with the builder. */
export function campaignCatalog(options: CampaignCatalogOptions): CampaignCatalog {
  const family = FAMILIES[options.family];
  // the last list read: the campaign-field columns and groupings are made of it, and stay while a refresh runs
  const known = signal<readonly AnyCampaign[]>([]);
  const rights = signal<ReadonlyMap<string, boolean>>(new Map());
  const access = new WeakMap<AnyCampaign, ReadonlySignal<CatalogAccess>>();
  let modified = new Map<string, number>();
  let freshRead = false;
  let reads = 0;
  const source = new AsyncValue(() => load(), {debounceMs: 0});
  // a list read with the same fields keeps the same column objects, and with them the table's rows
  const fields = computed(() => JSON.stringify(fieldsOf(known.value)));
  const table = new CatalogTable<AnyCampaign>({
    name: 'campaigns',
    label: 'Campaigns',
    items: source,
    keyOf: (c) => c.name,
    columns: computed(() => {
      fields.value;
      return catalogColumns(options.family, known.peek(), (c) => lastModified(c), options.narrow?.value);
    }),
    groupings: computed(() => {
      fields.value;
      return catalogGroupings(known.peek());
    }),
    defaultColumns: family.defaultColumns,
    defaultOrder: (a, b) => lastModified(b) - lastModified(a) || CODES.compare(a.name, b.name),
    emptyGroupLabel: i18n.noInformation,
    actionsVisible: 'always',
    access: (c) => accessOf(c),
    viewOnly: {icon: 'eye', tooltip: 'View only: you can open this campaign but not change it'},
    actions: (c) => actions(c),
    onActivate: (c) => options.open(c.name),
    empty: 'No campaigns yet. Start the first one under New Campaign.',
    footer: (shown) => {
      const viewOnly = shown.filter((c) => accessOf(c).value === 'view').length;
      return plural(shown.length, 'campaign', 'campaigns') + (viewOnly > 0 ? ` · ${viewOnly} view only` : '');
    },
    refresh: () => reload(true),
    persist: catalogPersistence(options.family),
  });
  return {
    table,
    chrome: () => table.chrome({columns: 'Edit visible columns in campaigns table',
      group: (g) => `Group Campaigns. Current: ${g}`, refresh: 'Refresh campaigns table'}),
    reload,
  };

  function reload(fresh = false): Promise<void> {
    freshRead ||= fresh;
    return source.refresh();
  }

  async function load(): Promise<readonly AnyCampaign[]> {
    const fresh = freshRead;
    freshRead = false;
    const generation = ++reads;
    const [campaigns, times] = await Promise.all([
      _package.loadCampaigns(options.app, options.deleted, fresh), readModified()]);
    const list = Object.values(campaigns) as AnyCampaign[];
    // a later read is under way: its answer is the one the table takes
    if (generation !== reads)
      return list;
    modified = times;
    known.value = list;
    void readRights(list, generation);
    options.onLoaded?.(list);
    return list;
  }

  /** When each campaign was saved last: its update_logs time (Hit Design, PeptiHit), else its campaign.json's; a
   * source that cannot be read leaves its campaigns to their Created date. */
  async function readModified(): Promise<Map<string, number>> {
    const app = options.app;
    const unread = (e: unknown) => {
      _package.logger.error(e);
      return null;
    };
    const [files, logs] = await Promise.all([_package.files.list(`${app}/campaigns`, true).catch(unread),
      family.saveLog ? api.queries.getLastModifiedAll(app).catch(unread) : null]);
    const times = new Map<string, number>();
    for (const file of files ?? []) {
      const time = file.updatedOn?.valueOf();
      if (file.name === CampaignJsonName && time !== undefined)
        times.set(file.path.split('/').slice(-2)[0], time);
    }
    if (logs) {
      const ids = logs.col('campaign_id')!;
      const at = logs.col('updated_at')!;
      for (let i = 0; i < logs.rowCount; i++) {
        if (!at.isNone(i))
          times.set(ids.get(i), at.get(i).valueOf());
      }
    }
    return times;
  }

  function lastModified(c: AnyCampaign): number {
    return modified.get(c.name) ?? (Date.parse(c.createDate) || 0);
  }

  /** The rights of every row land together, in place of the last answer for the rows it covered. */
  async function readRights(list: readonly AnyCampaign[], generation: number): Promise<void> {
    const me = DG.User.current().id;
    const membership: MembershipCache = new Map();
    const answers = await Promise.all(list.map(async (c) => [c.name,
      await checkEditPermissions(c.authorUserId ?? me, c.permissions ?? defaultPermissions, membership)] as const));
    if (generation === reads)
      rights.value = new Map(answers);
  }

  function accessOf(c: AnyCampaign): ReadonlySignal<CatalogAccess> {
    let a = access.get(c);
    if (a === undefined) {
      a = computed(() => {
        const edit = rights.value.get(c.name);
        return edit === undefined ? 'pending' : edit ? 'edit' : 'view';
      });
      access.set(c, a);
    }
    return a;
  }

  function actions(c: AnyCampaign): Action[] {
    return [
      {text: 'Manage campaign permissions', icon: 'users', requires: 'edit', run: () => void share(c)},
      {text: 'Delete campaign', icon: 'trash-alt', requires: 'edit', destructive: true,
        run: () => void remove(c)},
    ];
  }

  /** The landing's share writes campaign.json alone, under the campaign's lock and over what the server holds now
   * (another session may have saved since the list was read), and gives an authorless campaign its first author. */
  async function share(c: AnyCampaign): Promise<void> {
    await PermissionsEditor.show(c.permissions, {subject: subjectOf(c),
      author: c.authorUserFriendlyName, save: async (permissions) => {
        const app = options.app;
        const me = DG.User.current();
        const lock = await api.queries.acquireCampaignLock(app, c.name, me.friendlyName || me.name);
        if (!lock?.rowCount)
          throw new Error('Another session is saving this campaign. Try again in a moment.');
        try {
          const current = JSON.parse(await _package.files.readAsText(`${app}/campaigns/${c.name}/${CampaignJsonName}`));
          await _package.saveCampaignJson(app, {...current, permissions, authorUserId: current.authorUserId ?? me.id});
        } finally {
          await api.queries.releaseCampaignLock(app, c.name).catch((e) => _package.logger.error(e));
        }
        grok.shell.info('Permissions updated for campaign ' + options.label(c));
        await reload();
      }});
  }

  async function remove(c: AnyCampaign): Promise<void> {
    const deleted = await confirmDialog({title: 'Delete campaign', name: 'delete-campaign', confirmText: 'Delete',
      destructive: true, message: `Are you sure you want to delete campaign ${subjectOf(c, true)}?`,
      detail: `Its ${plural(c.rowCount ?? 0, 'molecule', 'molecules')} and its layout are removed for everyone. ` +
        'This cannot be undone.',
      onConfirm: () => options.remove(c.name)});
    if (!deleted)
      return;
    await reload();
    grok.shell.info(`Campaign ${options.label(c)} deleted`);
  }
}

/** The columns: Name first, then the family's getters, the optional Last Modified and every campaign field. */
export function catalogColumns(family: CampaignFamily, campaigns: readonly AnyCampaign[],
  modified: (c: AnyCampaign) => number, narrow = false): CatalogColumn<AnyCampaign>[] {
  const {getters, counts} = FAMILIES[family];
  const columns: CatalogColumn<AnyCampaign>[] = [{name: 'Name', value: getters['Name'], pinned: true,
    link: 'Continue Campaign', maxWidth: narrow ? '200px' : '320px'}];
  for (const name of Object.keys(getters).filter((g) => g !== 'Name')) {
    columns.push(counts.includes(name) ?
      {name, value: (c) => Number(getters[name](c)), align: 'right',
        description: name === 'Selected' ? 'Rows passing the saved filters' : undefined} :
      name === 'Status' ? {name, value: getters[name], render: (c) => statusCell(c.status)} :
        {name, value: getters[name]});
    if (name === LastModifiedByColumn) {
      const at = (c: AnyCampaign) => modified(c) > 0 ? new Date(modified(c)) : null;
      columns.push({name: LAST_MODIFIED, value: at, type: 'date', render: (c) => dayText(at(c))});
    }
  }
  for (const [field, templates] of fieldsOf(campaigns)) {
    columns.push({name: FIELD + field, header: field,
      value: (c) => fieldValue(c.campaignFields?.[field]),
      render: (c) => fieldText(c.campaignFields?.[field]), menuGroup: 'Campaign Fields',
      menuDetail: templates, maxWidth: '240px'});
  }
  return columns;
}

/** The groupings, named as `CampaignGrouping` and `campaignFields.<field>` (the kept choice), offered under the
 * captions of the columns they group by. A Template group is the template's name, trimmed and in any case: a
 * campaign without a template snapshot (LCNS-3) carries the name alone, where the others carry its key too. */
export function catalogGroupings(campaigns: readonly AnyCampaign[]): CatalogGrouping<AnyCampaign>[] {
  const template = (c: AnyCampaign) => (c.template?.name ?? c.templateName ?? '').trim();
  const groupings: CatalogGrouping<AnyCampaign>[] = [
    {name: CampaignGrouping.Template, key: (c) => template(c).toLowerCase(),
      groupLabel: (_key, first) => template(first)},
    {name: CampaignGrouping.Status, key: (c) => c.status},
    {name: CampaignGrouping.Author, key: (c) => c.authorUserFriendlyName},
    {name: CampaignGrouping.LastModifiedUser, label: LastModifiedByColumn, key: (c) => c.lastModifiedUserName},
  ];
  for (const [field, templates] of fieldsOf(campaigns)) {
    groupings.push({name: FIELD + field, label: field,
      key: (c) => fieldValue(c.campaignFields?.[field]), menuGroup: 'Campaign Fields',
      menuDetail: templates, groupLabel: (_key, first) => fieldText(first.campaignFields?.[field])});
  }
  return groupings;
}

/** The table's choices under the family's six keys, in their formats: the grouping as text ('None' for none), the
 * columns and the sort as JSON (no key for no sort), and the folded groups per grouping. The folds the platform
 * accordion kept per group (`Accordion:Hit-Design-campaigns-group-<key>`) are read while there are none. A save
 * writes only the keys whose choice changed, so a default nobody chose (the columns) is never frozen. */
export function catalogPersistence(family: CampaignFamily): CatalogPersistence {
  const {keys, defaultColumns} = FAMILIES[family];
  const json = (key: string): unknown => {
    try {
      return JSON.parse(localStorage.getItem(key) ?? 'null');
    }
    catch {
      return null;
    }
  };
  const texts = (s: CatalogState): [string, string | null][] => {
    const sort = SortOrder.parse(s.sort);
    return [
      [keys.grouping, s.groupBy ?? CampaignGrouping.None],
      [keys.columns, JSON.stringify(s.columns)],
      [keys.sorting, sort ? JSON.stringify({columnName: sort.column, ascending: !sort.descending}) : null],
      [keys.collapsed, JSON.stringify(s.collapsed)],
    ];
  };
  // what each key stands for now, as the table last had it
  const kept = new Map<string, string | null>();
  return {
    load: () => {
      const grouping = localStorage.getItem(keys.grouping);
      const sorting = json(keys.sorting) as {columnName?: unknown, ascending?: unknown} | null;
      const state: Partial<CatalogState> = {groupBy: grouping && grouping !== CampaignGrouping.None ? grouping : null,
        sort: typeof sorting?.columnName === 'string' && typeof sorting.ascending === 'boolean' ?
          SortOrder.of(sorting.columnName, !sorting.ascending) : ''};
      const columns = json(keys.columns);
      if (Array.isArray(columns) && columns.every((c) => typeof c === 'string'))
        state.columns = columns;
      state.collapsed = (json(keys.collapsed) ?? legacyCollapsed(state.groupBy ?? null)) as CatalogState['collapsed'];
      for (const [key, text] of texts({columns: defaultColumns, ...state} as CatalogState))
        kept.set(key, text);
      return state;
    },
    save: (s) => {
      for (const [key, text] of texts(s)) {
        if (kept.has(key) && kept.get(key) === text)
          continue;
        kept.set(key, text);
        if (text === null)
          localStorage.removeItem(key);
        else
          localStorage.setItem(key, text);
      }
    },
  };
}

/** A campaign field as the table sorts and groups it: a Date field, in whichever shape the forms wrote it
 * ({_value, isUtc, date} from the new-campaign form, {a, b, date} in older files), is a Date. An object of no known
 * shape is no value, never "[object Object]". */
export function fieldValue(value: unknown): unknown {
  if (value instanceof Date || value === null || typeof value !== 'object')
    return value;
  const raw = (value as {date?: unknown}).date ?? (value as {_value?: unknown})._value;
  const date = typeof raw === 'string' || typeof raw === 'number' ? new Date(raw) : null;
  return date !== null && !isNaN(date.getTime()) ? date : undefined;
}

/** A campaign field as its cell shows it: a date as the table's days are written (2026/3/30). */
export function fieldText(value: unknown): string {
  const v = fieldValue(value);
  return v instanceof Date ? dayText(v) : v === null || v === undefined ? '' : String(v);
}

/** A day as the campaigns write their Created date: 2026/3/30; no date, no text. */
export function dayText(d: Date | null): string {
  return d === null || isNaN(d.getTime()) ? '' : `${d.getFullYear()}/${d.getMonth() + 1}/${d.getDate()}`;
}

/** A status as a coloured dot and its text. */
function statusCell(status: string | undefined): HTMLElement | string {
  return status ? divH([dot(statusVariant(status)), span(status)], {gap: 's'}) : '';
}

function statusVariant(status: string): BadgeVariant {
  const s = status.toLowerCase();
  return s.includes('progress') ? 'accent' : s.includes('submitted') || s.includes('completed') ? 'success' :
    s.includes('hold') || s.includes('queued') || s.includes('awaiting') ? 'warning' : 'default';
}

/** Every field of the campaigns, by name, with the templates that define it. */
function fieldsOf(campaigns: readonly AnyCampaign[]): [string, string][] {
  const owners = new Map<string, Set<string>>();
  for (const c of campaigns) {
    for (const field of Object.keys(c.campaignFields ?? {}))
      owners.set(field, (owners.get(field) ?? new Set()).add(c.template?.name ?? c.templateName));
  }
  return [...owners.keys()].sort().map((field) => [field, [...owners.get(field)!].join(', ')]);
}

/** Groups start open now: of the accordion's folds, only a group closed on purpose stays closed. */
function legacyCollapsed(groupBy: string | null): CatalogState['collapsed'] {
  if (groupBy === null)
    return {};
  const folded = Object.keys(localStorage).filter((key) => key.startsWith(LEGACY_GROUP_KEY)).filter((key) => {
    try {
      return JSON.parse(localStorage.getItem(key)!)?.[key.slice(LEGACY_GROUP_KEY.length)] === false;
    }
    catch {
      return false;
    }
  });
  return {[groupBy]: folded.map((key) => key.slice(LEGACY_GROUP_KEY.length))};
}

/** "KHL-3 (BTK covalent follow-up)", or with `nameFirst` "BTK covalent follow-up (KHL-3)"; the code alone when it
 * is the name. */
function subjectOf(c: AnyCampaign, nameFirst = false): string {
  const name = c.friendlyName;
  if (!name || name === c.name)
    return c.name;
  return nameFirst ? `${name} (${c.name})` : `${c.name} (${name})`;
}
