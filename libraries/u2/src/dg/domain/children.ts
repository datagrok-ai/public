/* `domains.children` — the tables that refer to a source's rows, one tab each, over the current
   row: the child collection queried by the FK, New pre-filled with it, every child source in the
   parent's session — so a draft parent and its child drafts save as one transaction (the child's
   FK is the parent's draft id; the server orders the inserts). Once the parent is saved, the
   current row re-points to its real id and the tab re-queries by it. */
import {Control} from '../../core/component.js';
import {Scope} from '../../core/scope.js';
import {signal, untracked, ReadonlySignal} from '../../core/signals.js';
import {div, divH, span} from '../../core/elements.js';
import {TabStrip} from '../../components/containers/tabs.js';
import {loader} from '../../components/display/async-view.js';
import type {DomainSource} from '../../sources/domain-source.js';
import {Rows} from '../../sources/rows-like.js';
import {domains, DomainTable} from './index.js';
import {DomainApp} from './app.js';
import {DomainErrors} from './errors.js';

export type DomainChildrenMode = 'grid' | 'list';

export interface DomainChildrenOptions {
  /** The child tables to show (`<schema>.<table>` or the bare name); every one by default. */
  tables?: string[];
  /** `grid` (default): the platform grid; `list`: a list beside a form. */
  mode?: DomainChildrenMode;
}

/** One tab: the child table, the FK column pointing at the parent, the tab's id and label. */
export interface ChildEntry {
  id: string;
  table: DomainTable;
  fk: string;
  label: string;
}

export class DomainChildren extends Control {
  readonly tabs: TabStrip;
  readonly mode: DomainChildrenMode;
  /** The tabs, once the child tables are known — the parent loaded and the handles acquired. */
  readonly entries: ReadonlySignal<readonly ChildEntry[]>;

  private readonly _entries = signal<readonly ChildEntry[]>([]);
  private readonly _status: HTMLElement;
  private readonly _sources = new Map<string, DomainSource>();

  constructor(readonly parent: DomainSource, options: DomainChildrenOptions = {}) {
    super();
    this.mode = options.mode ?? 'grid';
    this.entries = this._entries;
    this.root.classList.add('u2-domain-children');
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
      void this._build(options.tables);
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
    let handles: DomainTable[];
    try {
      handles = await Promise.all(refs.map((c) => domains.table(`${c.schema}.${c.table}`)));
    } catch (e) {
      this._status.replaceChildren(span(DomainErrors.message(e), 'u2-domain-children-error'));
      return;
    }
    if (this.scope.isDisposed)
      return;
    this._status.replaceChildren();
    const entries = refs.map((c, i): ChildEntry => {
      const address = `${c.schema}.${c.table}`;
      // two FKs to one parent: two tabs for one table, told apart by the FK's caption
      const ambiguous = refs.filter((o) => `${o.schema}.${o.table}` === address).length > 1;
      const plural = DomainApp.titleOf(handles[i].info);
      return {id: `${address}.${c.fkColumn}`, table: handles[i], fk: c.fkColumn,
        label: ambiguous ? `${plural} (${c.label})` : plural};
    });
    for (const entry of entries) {
      const pane = div([], 'u2-domain-children-pane');
      pane.dataset.u2Part = 'pane';
      this.tabs.addTab({id: entry.id, label: entry.label, content: () => {
        this._watch(entry, pane);
        return pane;
      }});
    }
    this._entries.value = entries;
    void this._openFirstFilled(entries);
  }

  /** A row opens on a tab that has something in it: the first child table with rows, the first
   * tab when none has or the counts do not answer. A tab shown meanwhile is the user's. */
  private async _openFirstFilled(entries: readonly ChildEntry[]): Promise<void> {
    const row = this.parent.currentRow.peek();
    const opened = this.tabs.activeTab.peek();
    if (entries.length < 2 || row === null || Rows.isDraft(row))
      return;
    // a count nobody asked for: a table that refuses it just does not win the first tab
    const counts = await Promise.all(entries.map((e) =>
      e.table.table.count({filter: `${e.fk} = "${row.id}"`}).catch(() => 0)));
    const at = counts.findIndex((n) => n > 0);
    if (at > 0 && !this.scope.isDisposed && this.tabs.activeTab.peek() === opened)
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
      if (row === null) {
        const singular = this.parent.schema.info.singularName.toLowerCase() || 'row';
        pane.replaceChildren(span(`Select a ${singular}.`, 'u2-domain-children-hint'));
        return;
      }
      // everything built here — the source, the controls — is adopted by the pane's scope
      const pending = scope = new Scope();
      Scope.runWith(pending, () => untracked(() => {
        const draft = Rows.isDraft(row);
        const child = entry.table.source({query: draft ? '' : `${entry.fk} = "${row.id}"`,
          defaults: {[entry.fk]: row.id}, session: this.parent.session, empty: draft});
        this._sources.set(entry.id, child);
        // the parent's id is the same on every row here, and the column is not the user's to edit
        const body = this.mode === 'grid' ? domains.grid(child, {hiddenColumns: [entry.fk]}).root :
          divH([domains.list(child).root, domains.form(child).root], 'u2-domain-children-split');
        const empty = span(`No ${entry.label.toLowerCase()}.`, 'u2-domain-children-empty');
        pane.replaceChildren(divH([domains.newButton(child).root], 'u2-domain-children-toolbar'), body, empty);
        pending.effect(() => {
          const none = child.state.value === 'ready' && child.rows.items.value.length === 0;
          empty.hidden = !none;
          body.hidden = none;
        });
      }));
    });
  }
}
