import * as grok from 'datagrok-api/grok';
import {batch, computed, iconButton, signal, Dialog, Filters, Menu, TextInput, notify} from '@datagrok-libraries/u2';
import type {Signal, FilterProperty, FilterSchema} from '@datagrok-libraries/u2';
import type {DomainSource} from '../core/domain-source.js';
import {confirmDiscard} from '../core/session.js';

export interface DomainSearch {
  q: string;
  search: string;
  sort: string;
  at?: number;
}

const RECENT = 'ems-recent-searches';
const savedKey = (table: string) => `ems-searches|${table}`;
const same = (a: DomainSearch, b: DomainSearch) => a.q === b.q && a.search === b.search && a.sort === b.sort;
const label = (s: DomainSearch) => [s.q, s.search && `Search "${s.search}"`].filter((x) => x).join(' · ');
/** grok.userSettings keeps writes in memory until its 10 s sync, so a reload or a closed tab right after an
 * explicit save, delete or clear would lose it: send them now. Recent searches, written on every list change,
 * are left to the periodic sync. */
const persist = () => void grok.userSettings.flush();

/** The keys of `sort` over columns the table describes (a column the caller may not read, or one removed since, is
 * refused by the server for the whole read, so it is dropped and the rest of the order kept). */
export function sortable(sort: string, properties: readonly {name?: string}[]): string {
  const columns = new Set(properties.map((p) => p.name));
  return sort.split(',').filter((key) => key !== '' && columns.has(key.replace(/^!/, ''))).join(',');
}

/** The query text as the box reads it, in the form the box sends: a link's, a saved or a recent search's text goes
 * through the schema's grammar (case-insensitive and `[bracket]` names, its `like` escapes) and is written back
 * formatted, so the server parses what the box validated. Text the grammar or the schema refuses stays as it is: the
 * box shows its refusal and sends nothing. */
export function normalized(q: string, schema: {properties: readonly FilterProperty[]}): string {
  if (q.trim() === '' || schema.properties.length === 0)
    return q;
  const parsed = Filters.parse(q, schema as FilterSchema, 'domain');
  return parsed.problems.length > 0 ? q : Filters.format(parsed.root);
}

export class SavedSearches {
  static readonly changed = signal(0);
  private static readonly _current = new WeakMap<DomainSource, Signal<string | null>>();

  static of(source: DomainSource): DomainSearch {
    const q = source.query.value;
    return {q: typeof q === 'string' ? q : Filters.format(q), search: source.search.value, sort: source.sort.value};
  }

  static saved(table: string): Map<string, DomainSearch> {
    return new Map(Object.entries(grok.userSettings.get(savedKey(table)) ?? {}).map(([k, v]) => [k, JSON.parse(v)]));
  }

  /** The saved search `name` collides with (names compare trimmed and ignoring case), or null. */
  static taken(table: string, name: string): string | null {
    const key = name.trim().toLowerCase();
    return [...SavedSearches.saved(table).keys()].find((k) => k.trim().toLowerCase() === key) ?? null;
  }

  static save(table: string, name: string, search: DomainSearch): void {
    grok.userSettings.add(savedKey(table), name, JSON.stringify(search));
    persist();
    SavedSearches.changed.value++;
  }

  static remove(table: string, name: string): void {
    grok.userSettings.delete(savedKey(table), name);
    persist();
    SavedSearches.changed.value++;
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
    persist();
  }

  static current(source: DomainSource): Signal<string | null> {
    let current = SavedSearches._current.get(source);
    if (current === undefined)
      SavedSearches._current.set(source, current = signal<string | null>(null));
    return current;
  }

  /** The status bar's name for the list: the saved search it was set to, "(edited)" once it moved off it. */
  static title(source: DomainSource): string | null {
    SavedSearches.changed.value;
    const name = SavedSearches.current(source).value;
    const saved = name === null ? undefined : SavedSearches.saved(source.table).get(name);
    return saved === undefined ? null : same(saved, SavedSearches.of(source)) ? name : `${name} (edited)`;
  }

  static apply(source: DomainSource, search: DomainSearch): void {
    // a saved sort may name a column removed since: dropped, and said once, never written into the address
    const sort = sortable(search.sort, source.schema.properties);
    void confirmDiscard(source.session, {action: 'change the filter'}).then((ok) => {
      if (!ok)
        return;
      if (sort !== search.sort) {
        const gone = search.sort.split(',').filter((k) => k !== '' && !sort.split(',').includes(k))
          .map((k) => k.replace(/^!/, ''));
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
    const match = [...SavedSearches.saved(source.table)].find(([, s]) => same(s, search));
    if (match !== undefined)
      current.value = match[0];
    SavedSearches.record(source.table, search);
  }

  static historyButton(source: DomainSource): HTMLButtonElement {
    return iconButton('history', () => {}, {tooltip: 'Saved and recent searches', menu: () => {
      const menu = new Menu();
      const saved = SavedSearches.saved(source.table);
      for (const [name, search] of saved)
        menu.item(name, () => SavedSearches.apply(source, search));
      if (saved.size > 0)
        menu.separator();
      for (const r of SavedSearches.recent(source.table))
        menu.item(label(r), () => SavedSearches.apply(source, r), {detail: new Date(r.at!).toLocaleTimeString()});
      return menu.separator().item('Clear recent searches', () => SavedSearches.clearRecent(source.table));
    }});
  }

  static saveButton(source: DomainSource): HTMLButtonElement {
    const disabled = computed(() => label(SavedSearches.of(source)) === '' ? 'Nothing to save' : null);
    const el: HTMLButtonElement = iconButton('save', () => {
      const name = SavedSearches.current(source).peek();
      if (name === null)
        return SavedSearches._saveDialog(source);
      new Menu()
        .item(`Save changes to "${name}"`, () => SavedSearches._saveAs(source, name))
        .item('Save as new search…', () => SavedSearches._saveDialog(source))
        .show({anchor: el});
    }, {tooltip: 'Save search', disabled});
    return el;
  }

  private static _saveDialog(source: DomainSource): void {
    const name = new TextInput({label: 'Name', value: label(SavedSearches.of(source)).slice(0, 40)});
    // saving under a taken name would silently replace that search: refused until the name changes
    const clash = (value: string) => {
      const other = SavedSearches.taken(source.table, value);
      return other === null ? null : `A saved search named ${other} exists.`;
    };
    name.addValidator(clash);
    Dialog.create('Save search').add(name)
      .okDisabled(computed(() => name.value.trim() === '' ? 'Enter a name' : clash(name.value)), {footer: true})
      .onOK(() => SavedSearches._saveAs(source, name.peek().trim()), {text: 'SAVE', primary: true})
      .show();
  }

  private static _saveAs(source: DomainSource, name: string): void {
    SavedSearches.save(source.table, name, SavedSearches.of(source));
    SavedSearches.current(source).value = name;
  }
}
