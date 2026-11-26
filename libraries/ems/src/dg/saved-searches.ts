import * as grok from 'datagrok-api/grok';
import {batch, computed, untracked, confirmDialog, iconButton, signal, Dialog, Filters, Menu, TextInput, notify}
  from '@datagrok-libraries/u2';
import type {Signal, FilterSchema, IProperty} from '@datagrok-libraries/u2';
import type {DomainSource} from '../core/domain-source.js';
import {confirmDiscard} from '../core/session.js';
import {DomainTable} from './facade.js';
import {DomainErrors} from './errors.js';

export interface DomainSearch {
  q: string;
  search: string;
  sort: string;
  at?: number;
}

const RECENT = 'ems-recent-searches';
const savedKey = (table: string) => `ems-searches|${table}`;
const same = (a: DomainSearch, b: DomainSearch, schema: FilterSchema = {properties: []}) =>
  normalized(a.q, schema) === normalized(b.q, schema) && a.search === b.search && a.sort === b.sort;
const label = (s: DomainSearch) => [s.q, s.search && `Search "${s.search}"`].filter((x) => x).join(' · ');

/** The keys of `sort` over columns the list can be ordered by: the server refuses a read that names an unknown one. */
export function sortable(sort: string, properties: readonly IProperty[]): string {
  const columns = new Set(properties.filter(DomainTable.isSortable).map((p) => p.name));
  return sort.split(',').map((key) => key.trim())
    .filter((key) => key !== '' && columns.has(key.replace(/^!/, '').replace(/^prop\./, ''))).join(',');
}

/** The query text in the form the box sends, so the server parses what the box validated. */
export function normalized(q: string, schema: FilterSchema): string {
  if (q.trim() === '' || schema.properties.length === 0)
    return q;
  const parsed = Filters.parse(q, schema, 'domain');
  return parsed.problems.length > 0 ? q : Filters.format(parsed.root);
}

export class SavedSearches {
  // per table, the server's answer plus this tab's writes, written key by key: another tab may have changed others
  private static readonly _saved = new Map<string, Signal<Map<string, DomainSearch>>>();
  // per table, bumped by every read and write: a read answering after a later one is dropped
  private static readonly _turn = new Map<string, number>();
  private static readonly _current = new WeakMap<DomainSource, Signal<string | null>>();

  static of(source: DomainSource): DomainSearch {
    const q = source.query.value;
    return {q: typeof q === 'string' ? q : Filters.format(q), search: source.search.value, sort: source.sort.value};
  }

  static saved(table: string): Map<string, DomainSearch> {
    return SavedSearches._of(table).value;
  }

  /** Reads the saved searches of `table` from the server again: what another tab saved or deleted shows. */
  static async reload(table: string): Promise<void> {
    const turn = SavedSearches._next(table);
    const stored: Record<string, string> | null = await grok.dapi.userDataStorage.get(savedKey(table));
    if (turn === SavedSearches._turn.get(table))
      SavedSearches._of(table).value = new Map(Object.entries(stored ?? {}).map(([k, v]) => [k, JSON.parse(v)]));
  }

  private static _next(table: string): number {
    const turn = (SavedSearches._turn.get(table) ?? 0) + 1;
    SavedSearches._turn.set(table, turn);
    return turn;
  }

  private static _of(table: string): Signal<Map<string, DomainSearch>> {
    let saved = SavedSearches._saved.get(table);
    if (saved === undefined) {
      SavedSearches._saved.set(table, saved = signal(new Map()));
      SavedSearches.reload(table).catch((e) => DomainErrors.report(e));
    }
    return saved;
  }

  /** The saved search `name` collides with (names compare trimmed and ignoring case), or null. */
  static taken(table: string, name: string): string | null {
    const key = name.trim().toLowerCase();
    return [...SavedSearches.saved(table).keys()].find((k) => k.trim().toLowerCase() === key) ?? null;
  }

  static async save(table: string, name: string, search: DomainSearch): Promise<void> {
    await grok.dapi.userDataStorage.postValue(savedKey(table), name, JSON.stringify(search));
    SavedSearches._next(table);
    const saved = SavedSearches._of(table);
    saved.value = new Map(saved.peek()).set(name, search);
  }

  static async remove(table: string, name: string): Promise<void> {
    await grok.dapi.userDataStorage.remove(savedKey(table), name);
    SavedSearches._next(table);
    const saved = SavedSearches._of(table);
    saved.value = new Map([...saved.peek()].filter(([key]) => key !== name));
  }

  static recent(table: string): DomainSearch[] {
    return JSON.parse(grok.userSettings.get(RECENT)?.[table] ?? '[]');
  }

  static record(table: string, search: DomainSearch): void {
    const list = [{...search, at: Date.now()}, ...SavedSearches.recent(table).filter((r) => !same(r, search))];
    grok.userSettings.add(RECENT, table, JSON.stringify(list.slice(0, 20)));
  }

  static clearRecent(table: string): void {
    grok.userSettings.delete(RECENT, table);
    void grok.userSettings.flush();
  }

  static current(source: DomainSource): Signal<string | null> {
    let current = SavedSearches._current.get(source);
    if (current === undefined)
      SavedSearches._current.set(source, current = signal<string | null>(null));
    return current;
  }

  /** The status bar's name for the list: the saved search it was set to, "(edited)" once it moved off it. */
  static title(source: DomainSource): string | null {
    const name = SavedSearches.current(source).value;
    const saved = name === null ? undefined : SavedSearches.saved(source.table).get(name);
    return saved === undefined ? null : same(saved, SavedSearches.of(source), source.schema) ? name :
      `${name} (edited)`;
  }

  static apply(source: DomainSource, search: DomainSearch): void {
    const sort = sortable(search.sort, source.schema.properties);
    void confirmDiscard(source.session, {action: 'change the filter'}).then((ok) => {
      if (!ok)
        return;
      const gone = search.sort.split(',').map((k) => k.trim()).filter((k) => k !== '' && !sort.split(',').includes(k))
        .map((k) => k.replace(/^!/, ''));
      if (gone.length > 0) {
        notify.warning(`Sort by ${gone.join(', ')} is no longer available; showing ${sort === '' ?
          'the default order' : `the order by ${sort}`}`);
      }
      batch(() => {
        source.query.value = normalized(search.q, source.schema);
        source.search.value = search.search;
        source.sort.value = sort;
      });
    });
  }

  static query({q, search, sort}: DomainSearch): string {
    return `?${new URLSearchParams(Object.entries({q, search, sort}).filter(([, v]) => v !== ''))}`;
  }

  /** An effect over the list: names the saved search it holds, however it got there, and records it in history. */
  static track(source: DomainSource): void {
    const search = SavedSearches.of(source);
    if (source.deleted.value === 'only')
      return;
    const current = SavedSearches.current(source);
    if (search.q === '' && search.search === '') {
      current.value = null;
      return;
    }
    // the list's own search is what is recorded: a fresh read of the saved ones is not a new search
    const saved = untracked(() => SavedSearches.saved(source.table));
    const match = [...saved].find(([, s]) => same(s, search, source.schema));
    if (match !== undefined)
      current.value = match[0];
    SavedSearches.record(source.table, search);
  }

  /** Opens on this tab's copy at once and reads the server's again: what another tab saved shows the next time. */
  static historyButton(source: DomainSource): HTMLButtonElement {
    return iconButton('history', () => void SavedSearches._fresh(source.table),
      {tooltip: 'Saved and recent searches', menu: () => SavedSearches._history(source)});
  }

  // what this tab and the others saved, as the server has it now; this tab's copy where the read fails
  private static _fresh(table: string): Promise<void> {
    return SavedSearches.reload(table).catch((e) => void DomainErrors.report(e));
  }

  private static _history(source: DomainSource): Menu {
    const menu = new Menu();
    const saved = SavedSearches.saved(source.table);
    for (const [name, search] of saved)
      menu.item(name, () => SavedSearches.apply(source, search));
    if (saved.size > 0) {
      menu.group('Delete saved search');
      for (const name of saved.keys()) {
        menu.item(name, async () => await confirmDialog({title: 'Delete saved search', destructive: true,
          confirmText: 'Delete', message: `Delete the saved search "${name}"?`}) &&
          SavedSearches.remove(source.table, name));
      }
      menu.endGroup().separator();
    }
    for (const r of SavedSearches.recent(source.table)) {
      const at = new Date(r.at!);
      const detail = at.toDateString() === new Date().toDateString() ? at.toLocaleTimeString() :
        at.toLocaleString(undefined, {dateStyle: 'medium', timeStyle: 'short'});
      menu.item(label(r), () => SavedSearches.apply(source, r), {detail});
    }
    return menu.separator().item('Clear recent searches', () => SavedSearches.clearRecent(source.table));
  }

  static saveButton(source: DomainSource): HTMLButtonElement {
    const disabled = computed(() => label(SavedSearches.of(source)) === '' ? 'Nothing to save' : null);
    const el: HTMLButtonElement = iconButton('save', () => void SavedSearches._fresh(source.table).then(() => {
      const name = SavedSearches.current(source).peek();
      if (name === null || !SavedSearches.saved(source.table).has(name))
        return SavedSearches._saveDialog(source);
      new Menu()
        .item(`Save changes to "${name}"`, () => void SavedSearches._saveAs(source, name))
        .item('Save as new search…', () => SavedSearches._saveDialog(source))
        .show({anchor: el});
    }), {tooltip: 'Save search', disabled});
    return el;
  }

  private static _saveDialog(source: DomainSource): void {
    const dialog = Dialog.create('Save search');
    const name = dialog.build(() =>
      new TextInput({label: 'Name', value: label(SavedSearches.of(source)).slice(0, 40)}));
    const clash = (value: string) => {
      const other = SavedSearches.taken(source.table, value);
      return other === null ? null : `A saved search named ${other} exists.`;
    };
    name.addValidator(clash);
    dialog.add(name)
      .okDisabled(computed(() => name.value.trim() === '' ? 'Enter a name' : clash(name.value)))
      .onOK(() => SavedSearches._saveAs(source, name.peek().trim()), {text: 'SAVE', primary: true})
      .show();
  }

  private static async _saveAs(source: DomainSource, name: string): Promise<void> {
    await SavedSearches.save(source.table, name, SavedSearches.of(source));
    SavedSearches.current(source).value = name;
  }
}
