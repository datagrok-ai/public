import * as grok from 'datagrok-api/grok';
import {batch, computed, iconButton, signal, Dialog, Filters, Menu, TextInput} from '@datagrok-libraries/u2';
import type {Signal} from '@datagrok-libraries/u2';
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

  static save(table: string, name: string, search: DomainSearch): void {
    grok.userSettings.add(savedKey(table), name, JSON.stringify(search));
    SavedSearches.changed.value++;
  }

  static remove(table: string, name: string): void {
    grok.userSettings.delete(savedKey(table), name);
    SavedSearches.changed.value++;
  }

  static recent(table: string): DomainSearch[] {
    return JSON.parse(grok.userSettings.get(RECENT)?.[table] ?? '[]');
  }

  static record(table: string, search: DomainSearch): void {
    const list = [{...search, at: Date.now()}, ...SavedSearches.recent(table).filter((r) => !same(r, search))];
    grok.userSettings.add(RECENT, table, JSON.stringify(list.slice(0, 20)));
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
    void confirmDiscard(source.session, {action: 'change the filter'}).then((ok) => ok && batch(() => {
      source.query.value = search.q;
      source.search.value = search.search;
      source.sort.value = search.sort;
    }));
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
      return menu.separator().item('Clear recent searches', () => grok.userSettings.delete(RECENT, source.table));
    }});
  }

  static saveButton(source: DomainSource): HTMLButtonElement {
    const disabled = computed(() => label(SavedSearches.of(source)) === '' ? 'Nothing to save' : false);
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
    Dialog.create('Save search').add(name)
      .okDisabled(computed(() => name.value.value.trim() === '' ? 'Enter a name' : null))
      .onOK(() => SavedSearches._saveAs(source, name.value.peek().trim()), {text: 'SAVE', primary: true})
      .show();
  }

  private static _saveAs(source: DomainSource, name: string): void {
    SavedSearches.save(source.table, name, SavedSearches.of(source));
    SavedSearches.current(source).value = name;
  }
}
