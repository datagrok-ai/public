import * as grok from 'datagrok-api/grok';
import type * as DG from 'datagrok-api/dg';
import {domains, DomainTable} from './facade.js';
import type {DomainAppViewOptions} from './facade.js';
import {DomainApp} from './app.js';
import {DomainAddress} from './address.js';
import {SavedSearches} from './saved-searches.js';
import {rawEffect} from 'datagrok-api/u2core';

/** What `DomainTable.app()` takes for one view, less what the views decide: its name and path. */
export type DomainViewOptions = Omit<DomainAppViewOptions, 'name' | 'path'>;

/** One table's view: its name — the tab's and the tree item's — and its segment under the app's path, the name
 * without spaces by default; '' is the app's path itself, where an app's main table keeps the URL it had. */
export interface DomainViewEntry {
  /** `<schema>.<table>`. */
  table: string;
  name: string;
  segment?: string;
  options?: DomainViewOptions;
}

export interface DomainViewsOptions {
  /** The app's route — what its app func is routed at, e.g. `/apps/Grit/Issues`. */
  path: string;
  /** The first is the view an address naming none opens. */
  views: DomainViewEntry[];
  /** The handle a view is made over — the app's configured one, whose validators and renderer then hold in it;
   * `domains.table` by default. */
  handle?: (address: string) => Promise<DomainTable>;
  /** Runs on each view's app once it is made — what the app wires onto it beyond the options. */
  ready?: (app: DomainApp, table: string) => void;
}

/** Several tables as one app: a view per table at `<path>/<Segment>`, and a row of one of the tables opens from
 * anywhere in the app's own view of the table, not the platform's entity view. */
export class DomainViews {
  // where each table's rows open: the views made over it last
  private static readonly _byTable = new Map<string, DomainViews>();
  private static _deleteSaved: {unsubscribe(): void} | undefined;

  constructor(private readonly _options: DomainViewsOptions) {
    for (const entry of _options.views)
      DomainViews._byTable.set(entry.table, this);
  }

  /** Opens a saved row of `table` in the view of it some views have — brought forward at the row, or docked; false
   * where no views show the table. */
  static openRow(table: string, id: string): boolean {
    const views = DomainViews._byTable.get(table);
    if (views !== undefined)
      void views.show(`${views._segment(views._entry(table))}?entity=${encodeURIComponent(DomainAddress.escape(id))}`);
    return views !== undefined;
  }

  /** What the app func answers the router for an address under the path — `Batches`, `Compounds?entity=<id>`; '' is
   * the first view's list: a new view, which the router docks, or null where the open view of the table takes the
   * address and comes forward. Back and Forward inside an open view reach the app func too (the router cannot ask a
   * plain view for its path), and the open view answered again would be docked twice. */
  async open(address = ''): Promise<DG.ViewBase | null> {
    const {entry, rest} = this._route(address);
    return DomainApp.openAt(this.baseOf(entry.table), this._at(entry, rest)) === undefined ? this._new(address) : null;
  }

  /** An address from inside the app — a tree item, a link: the open view comes forward at it, else a new one
   * docks. */
  async show(address: string): Promise<DG.ViewBase> {
    const {entry, rest} = this._route(address);
    return DomainApp.openAt(this.baseOf(entry.table), rest) ?? grok.shell.addView(await this._new(address));
  }

  /** The create page of a table's view — the open view's, else a new one's. */
  async create(table: string): Promise<void> {
    await DomainApp.of(await this.show(this._segment(this._entry(table))))!.create();
  }

  /** A new view of a table's list, not docked. */
  async view(table: string): Promise<DG.ViewBase> {
    const entry = this._entry(table);
    const handle = await (this._options.handle ?? domains.table)(table);
    const view = handle.app({...entry.options, name: entry.name, path: this.baseOf(table)});
    this._options.ready?.(DomainApp.of(view)!, table);
    return view;
  }

  static async reload(schema: string): Promise<void> {
    DomainTable.forget(schema);
    for (const app of [...DomainApp.live]) {
      const old = DomainApp.viewOf(app);
      if (!app.table.address.startsWith(`${schema}.`) || old === undefined || app.session.isDirty.peek())
        continue;
      const views = DomainViews._byTable.get(app.table.address);
      const view = views !== undefined ? await views._new(app.path.peek().slice(views._options.path.length)) :
        (await domains.table(app.table.address)).app({path: app.base});
      grok.shell.addView(view);
      old.close();
    }
  }

  /** An item per view under `node` — the app's tree browser: one selected brings its view forward at the list —
   * then the user's saved searches. */
  tree(node: DG.TreeViewGroup): void {
    for (const entry of this._options.views)
      node.item(entry.name).onSelected.subscribe(() => void this.show(this._segment(entry)));
    const saved = node.group('Saved searches', null, false);
    rawEffect(() => {
      SavedSearches.changed.value;
      for (const child of [...saved.children])
        child.remove();
      for (const entry of this._options.views) {
        for (const [name, search] of SavedSearches.saved(entry.table)) {
          saved.item(`${entry.name}: ${name}`, {savedSearch: {table: entry.table, name}})
            .onSelected.subscribe(() => void this.show(`${this._segment(entry)}${SavedSearches.query(search)}`));
        }
      }
    });
    type ContextMenu = {args: {item?: DG.TreeViewNode, menu: DG.Menu}};
    DomainViews._deleteSaved ??= grok.events.onContextMenu.subscribe((e: ContextMenu) => {
      const saved = e.args.item?.value?.savedSearch;
      if (saved !== undefined)
        e.args.menu.item('Delete', () => SavedSearches.remove(saved.table, saved.name));
    });
  }

  /** The route of a table's view: the app's path and the view's segment. */
  baseOf(table: string): string {
    const segment = this._segment(this._entry(table));
    return segment === '' ? this._options.path : `${this._options.path}/${segment}`;
  }

  private async _new(address: string): Promise<DG.ViewBase> {
    const {entry, rest} = this._route(address);
    const view = await this.view(entry.table);
    void DomainApp.of(view)!.open(this._at(entry, rest));
    return view;
  }

  // the router hands over the path alone: where the address bar is this view's, its query is replayed from there
  private _at(entry: DomainViewEntry, rest: string): string {
    const own = window.location.pathname.toLowerCase().endsWith(this.baseOf(entry.table).toLowerCase());
    return rest || (own ? `${window.location.pathname}${window.location.search}` : '');
  }

  private _entry(table: string): DomainViewEntry {
    return this._options.views.find((v) => v.table === table)!;
  }

  private _segment(entry: DomainViewEntry): string {
    return entry.segment ?? entry.name.replace(/ /g, '');
  }

  private _route(address: string): {entry: DomainViewEntry, rest: string} {
    const [, segment, rest] = /^\/?([^/?]*)(.*)$/.exec(address)!;
    const views = this._options.views;
    const entry = views.find((v) => this._segment(v).toLowerCase() === segment.toLowerCase());
    return entry === undefined ? {entry: views[0], rest: address.startsWith('?') ? address : ''} : {entry, rest};
  }
}
