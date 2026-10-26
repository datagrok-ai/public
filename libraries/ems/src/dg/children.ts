/* `domains.children` — the tables that refer to a source's rows, one tab each, over the current
   row: the child collection queried by the FK, New pre-filled with it, every child source in the
   parent's session — so a draft parent and its child drafts save as one transaction (the child's
   FK is the parent's draft id; the server orders the inserts). Once the parent is saved, the
   current row re-points to its real id and the tab re-queries by it. */
import {Control, Scope, signal, computed, untracked, ReadonlySignal, Signal, div, divH, span, TabStrip, loader,
  Rows} from '@datagrok-libraries/u2';
import type {RowView} from '@datagrok-libraries/u2';
import {DomainSource} from '../core/domain-source.js';
import {toPlural} from '../core/domain-backend.js';
import {domains, DomainTable} from './facade.js';
import {DomainApp} from './app.js';
import {DomainErrors} from './errors.js';

export type DomainChildrenMode = 'grid' | 'list';

export interface DomainChildrenOptions {
  /** The child tables to show (`<schema>.<table>` or the bare name); every one by default. */
  tables?: string[];
  /** `grid` (default): the platform grid; `list`: a list beside a form. */
  mode?: DomainChildrenMode;
  /** Per child table (`<schema>.<table>` or the bare name), in grid mode: the columns the grid shows,
   * in order, and the column whose cell opens the child's own page. */
  grids?: Record<string, {columns?: string[], open?: string}>;
}

/** One tab: the child table, the FK column pointing at the parent, the tab's id and label. */
export interface ChildEntry {
  id: string;
  table: DomainTable;
  fk: string;
  label: string;
}

/** What a tab follows: its rows' count, the child source while the tab is shown, and whether the
 * pane holds nothing but its empty message. */
interface TabState {
  count: Signal<number | null>;
  current: Signal<DomainSource | undefined>;
  blank: Signal<boolean>;
}

export class DomainChildren extends Control {
  readonly tabs: TabStrip;
  readonly mode: DomainChildrenMode;
  /** The tabs, once the child tables are known — the parent loaded and the handles acquired; a
   * table that did not open has a tab with its error, and no entry. */
  readonly entries: ReadonlySignal<readonly ChildEntry[]>;

  private readonly _entries = signal<readonly ChildEntry[]>([]);
  private readonly _status: HTMLElement;
  private readonly _sources = new Map<string, DomainSource>();
  private readonly _tabStates = new Map<string, TabState>();
  private readonly _grids: NonNullable<DomainChildrenOptions['grids']>;

  constructor(readonly parent: DomainSource, options: DomainChildrenOptions = {}) {
    super();
    this.mode = options.mode ?? 'grid';
    this._grids = options.grids ?? {};
    this.entries = this._entries;
    this.root.classList.add('u2-domain-children', `u2-domain-children-${this.mode}`);
    this.root.dataset.u2 = 'domain-children';
    this._status = div([], 'u2-domain-children-status');
    this._status.dataset.u2Part = 'status';
    this.tabs = this.runInScope(() => new TabStrip());
    this.tabs.root.classList.add('u2-domain-children-tabs');
    this.root.append(this._status, this.tabs.root);
    // the child tables are known with the parent's table: build the tabs once it is loaded
    let built = false;
    this.effect(() => {
      if (built || parent.state.value !== 'ready')
        return;
      built = true;
      // built once and kept: owned by this control, not by the effect run that triggered it
      void this.runInScope(() => this._build(options.tables));
    });
  }

  /** The source behind a built tab (its id is `entries[i].id`); undefined until the tab was shown. */
  child(id: string): DomainSource | undefined {
    return this._sources.get(id);
  }

  private async _build(only?: string[]): Promise<void> {
    const at = (c: {schema: string, table: string}) =>
      Math.max(only!.indexOf(`${c.schema}.${c.table}`), only!.indexOf(c.table));
    const refs = this.parent.schema.info.childTables.filter((c) =>
      only === undefined || only.includes(`${c.schema}.${c.table}`) || only.includes(c.table));
    // the tabs are the caller's order where it named them; the registry's otherwise
    if (only !== undefined)
      refs.sort((a, b) => at(a) - at(b));
    this._status.replaceChildren(loader());
    const opened = await Promise.allSettled(refs.map((c) => domains.table(`${c.schema}.${c.table}`)));
    if (this.scope.isDisposed)
      return;
    this._status.replaceChildren();
    const entries: ChildEntry[] = [];
    for (const [i, c] of refs.entries()) {
      const address = `${c.schema}.${c.table}`;
      const id = `${address}.${c.fkColumn}`;
      const handle = opened[i];
      const title = DomainApp.titleOf(handle.status === 'fulfilled' ? handle.value.info :
        {pluralName: toPlural(c.table)});
      // two FKs to one parent: two tabs for one table, told apart by the FK's caption
      const ambiguous = refs.filter((o) => `${o.schema}.${o.table}` === address).length > 1;
      const label = ambiguous ? `${title} (${c.label})` : title;
      const pane = div([], 'u2-domain-children-pane');
      pane.dataset.u2Part = 'pane';
      // a child table that does not open fails its own tab; the others still work
      if (handle.status === 'rejected') {
        pane.append(span(DomainErrors.message(handle.reason), 'u2-domain-children-error'));
        this.tabs.addTab({id, label, content: () => pane});
        continue;
      }
      const entry: ChildEntry = {id, table: handle.value, fk: c.fkColumn, label};
      entries.push(entry);
      const state: TabState = {count: signal(null), current: signal(undefined), blank: signal(false)};
      this._tabStates.set(id, state);
      const validity = computed(() => DomainChildren._problem(state.current.value));
      this.tabs.addTab({id, label, validity, content: () => {
        this._watch(entry, pane);
        return pane;
      }});
      const counted = computed(() => state.count.value === null ? label : `${label} (${state.count.value})`);
      this.effect(() => {
        const text = counted.value;
        untracked(() => this.tabs.setTabLabel(id, text));
      });
      this.effect(() => {
        const dirty = state.current.value?.isDirty.value ?? false;
        untracked(() => this.tabs.setDirty(id, dirty));
      });
    }
    this.effect(() => {
      const shown = this.tabs.activeTab.value;
      const blank = shown !== null && this._tabStates.get(shown)?.blank.value === true;
      this.root.classList.toggle('u2-domain-children-blank', blank);
    });
    this._entries.value = entries;
    let picked = false;
    let countedFor: string | null | undefined;
    this.effect(() => {
      const row = this.parent.currentRow.value;
      // a saved draft re-points the row mid-save: the counts follow once it has landed
      if ((row?.id ?? null) === countedFor || this.parent.session.isSaving.value)
        return;
      countedFor = row?.id ?? null;
      untracked(() => void this._count(entries, row, !picked));
      picked = true;
    });
  }

  /** What a tab is marked invalid with: the refusal its child source carries, never a failed load,
   * which the pane shows itself. */
  private static _problem(source: DomainSource | undefined): string | null {
    const error = source?.error.value;
    const code = DomainErrors.codeOf(error);
    return error === undefined || source!.state.value === 'error' ||
      (code !== 'validation' && code !== DomainSource.REFUSED) ? null : DomainErrors.message(error, source);
  }

  /** Every tab's rows counted for the parent row — its label says how many. A row opens on a tab
   * that has something in it (`pick`): the first child table with rows, the first tab when none has
   * or the counts do not answer. A tab shown meanwhile is the user's. */
  private async _count(entries: readonly ChildEntry[], row: RowView | null, pick: boolean): Promise<void> {
    for (const e of entries)
      this._tabStates.get(e.id)!.count.value = null;
    if (row === null || Rows.isDraft(row))
      return;
    const opened = this.tabs.activeTab.peek();
    // a count nobody asked for: a table that refuses it just goes without, and does not win the first tab
    const counts = await Promise.all(entries.map((e) =>
      e.table.table.count({filter: `${e.fk} = "${row.id}"`}).catch(() => null)));
    if (this.scope.isDisposed || this.parent.currentRow.peek()?.id !== row.id)
      return;
    for (const [i, e] of entries.entries())
      this._tabStates.get(e.id)!.count.value ??= counts[i];
    // the tab already shown may be one whose table did not open: it has no entry and no rows
    if (!pick || entries.every((e) => e.id === opened))
      return;
    const at = counts.findIndex((n) => n !== null && n > 0);
    if (at >= 0 && entries[at].id !== opened && this.tabs.activeTab.peek() === opened)
      this.tabs.activeTab.value = entries[at].id;
  }

  /** A shown tab follows the parent's current row: a new row is a new child source, the old one
   * disposed with everything built over it. */
  private _watch(entry: ChildEntry, pane: HTMLElement): void {
    let scope: Scope | undefined;
    let last: string | undefined;
    this.own(() => scope?.dispose());
    this.effect(() => {
      const row = this.parent.currentRow.value;
      // a saved draft re-points the row mid-save: the child (part of that batch) rebuilds after it
      if (row?.id === last || this.parent.session.isSaving.value)
        return;
      last = row?.id;
      scope?.dispose();
      scope = undefined;
      this._sources.delete(entry.id);
      const state = this._tabStates.get(entry.id)!;
      state.current.value = undefined;
      if (row === null) {
        pane.replaceChildren(span(`Select a ${DomainTable.noun(this.parent.schema.info, true)}.`,
          'u2-domain-children-hint'));
        state.blank.value = true;
        return;
      }
      // everything built here — the source, the controls — is adopted by the pane's scope
      const pending = scope = new Scope();
      Scope.runWith(pending, () => untracked(() => {
        const draft = Rows.isDraft(row);
        const child = entry.table.source({query: draft ? '' : `${entry.fk} = "${row.id}"`,
          defaults: {[entry.fk]: row.id}, session: this.parent.session, empty: draft});
        this._sources.set(entry.id, child);
        state.current.value = child;
        // the parent's id is the same on every row here, and the column is not the user's to edit
        const shape = this._grids[entry.table.address] ?? this._grids[entry.table.address.split('.').pop()!] ?? {};
        const open = shape.open;
        const grid = this.mode === 'grid' ? domains.grid(child, {hiddenColumns: [entry.fk], columns: shape.columns,
          ...(open === undefined ? {} : {open: {column: open, run: (r: RowView) => entry.table.open(r)}})}) : null;
        const body = grid?.root ??
          divH([domains.list(child).root, domains.form(child).root], 'u2-domain-children-split');
        const empty = span(`No ${entry.label.toLowerCase()}.`, 'u2-domain-children-empty');
        const add = domains.newButton(child);
        // the press leaves the focus where the new row is typed: the grid's row; a list's form takes it itself
        if (grid !== null)
          add.root.addEventListener('click', () => queueMicrotask(() => grid.grid.root.focus()));
        pane.replaceChildren(divH([add.root], 'u2-domain-children-toolbar'), body, empty);
        pending.effect(() => {
          const ready = child.state.value === 'ready';
          const none = ready && child.rows.items.value.length === 0;
          empty.hidden = !none;
          body.hidden = none;
          state.blank.value = none;
          if (ready && !draft)
            state.count.value = child.total.value;
        });
      }));
    });
  }
}
