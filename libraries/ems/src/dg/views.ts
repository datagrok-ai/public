/* Several tables as one app: a view per table at `<path>/<Segment>` — what the app func answers the router, a
   deep link (`?entity=<id>`, `?q=`) replayed by the view — an item per view in the app's tree, and where a row of
   one of the tables opens from anywhere: the app's own view of the table, not the platform's entity view. */
import * as grok from 'datagrok-api/grok';
import type * as DG from 'datagrok-api/dg';
import {domains} from './facade.js';
import type {DomainAppViewOptions, DomainTable} from './facade.js';
import {DomainApp} from './app.js';
import {DomainAddress} from './address.js';

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

export class DomainViews {
  /** Where each table's rows open: every views object over it, the latest last. */
  private static readonly _byTable = new Map<string, DomainViews[]>();

  constructor(private readonly _options: DomainViewsOptions) {
    for (const entry of _options.views)
      DomainViews._byTable.set(entry.table, [...DomainViews._byTable.get(entry.table) ?? [], this]);
  }

  /** Opens a saved row of `table` in the view of it some views have — brought forward at the row, or docked; false
   * where no views show the table. Views with an app open over the table win over the latest made. */
  static openRow(table: string, id: string): boolean {
    const all = DomainViews._byTable.get(table) ?? [];
    const open = all.filter((v) => [...DomainApp.live].some((app) => app.base === v.baseOf(table)));
    const views = open[open.length - 1] ?? all[all.length - 1];
    if (views !== undefined)
      void views.show(`${views._segment(views._entry(table))}?entity=${encodeURIComponent(DomainAddress.escape(id))}`);
    return views !== undefined;
  }

  /** Takes these views out of where rows open. */
  dispose(): void {
    for (const entry of this._options.views) {
      const left = (DomainViews._byTable.get(entry.table) ?? []).filter((v) => v !== this);
      if (left.length === 0)
        DomainViews._byTable.delete(entry.table);
      else
        DomainViews._byTable.set(entry.table, left);
    }
  }

  /** A new view at an address under the path — `Batches`, `Compounds?entity=<id>`; '' is the first view's list —
   * what the app func answers the router, which docks it. The router hands over the path alone: where the address
   * bar is this view's, its query is replayed by the view itself. */
  async open(address = ''): Promise<DG.ViewBase> {
    const {entry, rest} = this._route(address);
    const view = await this.view(entry.table);
    const own = window.location.pathname.toLowerCase().endsWith(this.baseOf(entry.table).toLowerCase());
    void DomainApp.of(view)!.open(rest || (own ? undefined : ''));
    return view;
  }

  /** An address from inside the app — a tree item, a link: the open view comes forward at it, else a new one
   * docks. */
  async show(address: string): Promise<DG.ViewBase> {
    const {entry, rest} = this._route(address);
    return DomainApp.openAt(this.baseOf(entry.table), rest) ?? grok.shell.addView(await this.open(address));
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

  /** An item per view under `node` — the app's tree browser: one selected brings its view forward at the list. */
  tree(node: DG.TreeViewGroup): void {
    for (const entry of this._options.views)
      node.item(entry.name).onSelected.subscribe(() => void this.show(this._segment(entry)));
  }

  /** The route of a table's view: the app's path and the view's segment. */
  baseOf(table: string): string {
    const segment = this._segment(this._entry(table));
    return segment === '' ? this._options.path : `${this._options.path}/${segment}`;
  }

  private _entry(table: string): DomainViewEntry {
    return this._options.views.find((v) => v.table === table)!;
  }

  private _segment(entry: DomainViewEntry): string {
    return entry.segment ?? entry.name.replace(/ /g, '');
  }

  /** The view an address names by its first segment — the first view where it names none — and the rest of it. */
  private _route(address: string): {entry: DomainViewEntry, rest: string} {
    const [, segment, rest] = /^\/?([^/?]*)(.*)$/.exec(address)!;
    const views = this._options.views;
    const entry = views.find((v) => this._segment(v).toLowerCase() === segment.toLowerCase());
    return entry === undefined ? {entry: views[0], rest: address.startsWith('?') ? address : ''} : {entry, rest};
  }
}
