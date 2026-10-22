/* Several tables as one app: a view per table at `<path>/<Segment>` — what the app func answers the router, a
   deep link (`?entity=<id>`, `?entity=new`, `?q=`) replayed by the view — an item per view in the app's tree, and
   where a row of one of the tables opens when no view shows it: a reference's open control and a list's Open land
   in the app's own view of the table, not the platform's entity view. */
import * as grok from 'datagrok-api/grok';
import type * as DG from 'datagrok-api/dg';
import {domains} from './facade.js';
import type {DomainAppViewOptions, DomainTable} from './facade.js';
import {DomainApp} from './app.js';
import {DomainAddress} from './address.js';
import type {ChildEntry} from './children.js';

/** What `DomainTable.app()` takes for one view, less what the views decide: its name and path. */
export type DomainViewOptions = Omit<DomainAppViewOptions, 'name' | 'path'>;

/** One table's view: its name — the tab's and the tree item's — and its segment under the app's path, the name
 * without spaces by default. With `childPages`, New on a saved row's child tab opens the child table's own view at a
 * create page, so every child table needs a view too. `options` are the view's, or built from the table's handle. */
export interface DomainViewEntry {
  /** `<schema>.<table>`. */
  table: string;
  name: string;
  segment?: string;
  childPages?: boolean;
  options?: DomainViewOptions | ((table: DomainTable) => DomainViewOptions | Promise<DomainViewOptions>);
}

export interface DomainViewsOptions {
  /** The app's route — what its app func is routed at, e.g. `/apps/Grit/Issues`. */
  path: string;
  /** The first is the view an address naming none opens. */
  views: DomainViewEntry[];
  /** The handle a view is made over — the app's configured one, whose derivations, validators and renderer then
   * hold in it; `domains.table` by default. */
  handle?: (address: string) => Promise<DomainTable>;
}

export class DomainViews {
  /** Where each table's rows open: the views made over it last. */
  private static readonly _byTable = new Map<string, DomainViews>();

  constructor(private readonly _options: DomainViewsOptions) {
    for (const entry of _options.views)
      DomainViews._byTable.set(entry.table, this);
  }

  /** Opens a saved row of `table` where an app shows it: the app open over the table, else the view of it the
   * views made over it have — brought forward at the row, or docked; false where neither is. */
  static openRow(table: string, id: string): boolean {
    if (DomainApp.activate(DomainApp.baseOf(table), id))
      return true;
    const views = DomainViews._byTable.get(table);
    if (views !== undefined)
      void views.show(`${views._segment(views._entry(table))}?entity=${encodeURIComponent(DomainAddress.escape(id))}`);
    return views !== undefined;
  }

  /** Whether {@link openRow} has an app to open a row of `table` in. */
  static opens(table: string): boolean {
    return DomainApp.baseOf(table) !== undefined || DomainViews._byTable.has(table);
  }

  /** A new view at an address under the path — `Batches`, `Compounds?entity=new`; '' is the first view's list —
   * what the app func answers the router, which docks it. The router hands over the path alone: a cold deep link's
   * query is replayed from the address bar by the view itself. */
  async open(address = ''): Promise<DG.ViewBase> {
    const {entry, rest} = this._route(address);
    const view = await this.view(entry.table);
    void DomainApp.of(view)!.open(rest);
    return view;
  }

  /** An address from inside the app — a tree item, a link: the open view comes forward at it, else a new one
   * docks. */
  async show(address: string): Promise<DG.ViewBase> {
    const {entry, rest} = this._route(address);
    return DomainApp.openAt(this.baseOf(entry.table), rest) ?? grok.shell.addView(await this.open(address));
  }

  /** The create page of a table's view over `values` — the open view's, else a new one's. */
  async create(table: string, values: Record<string, unknown> = {}): Promise<void> {
    if (!DomainApp.createAt(this.baseOf(table), values))
      await DomainApp.of(await this.show(this._segment(this._entry(table))))!.create(values);
  }

  /** A new view of a table's list, not docked. */
  async view(table: string): Promise<DG.ViewBase> {
    const entry = this._entry(table);
    const handle = await (this._options.handle ?? domains.table)(table);
    const options = typeof entry.options === 'function' ? await entry.options(handle) : entry.options ?? {};
    const children = !entry.childPages || options.children === false ? options.children :
      {...(options.children === true ? {} : options.children),
        onNew: (child: ChildEntry, values: Record<string, unknown>) => void this.create(child.table.address, values)};
    return handle.app({...options, children, name: entry.name, path: this.baseOf(table)});
  }

  /** An item per view under `node` — the app's tree browser: one selected brings its view forward at the list. */
  tree(node: DG.TreeViewGroup): void {
    for (const entry of this._options.views)
      node.item(entry.name).onSelected.subscribe(() => void this.show(this._segment(entry)));
  }

  /** The route of a table's view: the app's path and the view's segment. */
  baseOf(table: string): string {
    return `${this._options.path}/${this._segment(this._entry(table))}`;
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
    return {entry: views.find((v) => this._segment(v).toLowerCase() === segment.toLowerCase()) ?? views[0], rest};
  }
}
