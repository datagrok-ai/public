/* `domains.table(address)` — one await, then a table handle everything else hangs off (GOAL
   "What it looks like"): sources over the table, and the per-table registries an app fills once —
   actions (ribbon, row hover, context menu), the renderer (lists, pickers, chips) and validators
   (forms). Registers nothing with the platform: rendering reuses the handler `ObjectHandler.forEntity`
   already resolves for the table's rows. */
import * as DG from 'datagrok-api/dg';
import {Access} from '../../core/access.js';
import type {Capability} from '../../core/access.js';
import type {IProperty} from '../../core/property-like.js';
import type {ObjectRenderer} from '../../core/object-renderer.js';
import type {Action} from '../../components/actions/actions.js';
import {backends} from '../../sources/backends.js';
import type {DomainTableInfoLike, DomainTableLike} from '../../sources/domain-backend.js';
import {DomainSource} from '../../sources/domain-source.js';
import type {DomainSourceOptions} from '../../sources/domain-source.js';
import {Rows} from '../../sources/rows-like.js';
import type {RowView} from '../../sources/rows-like.js';
import type {ReadonlySignal} from '../../core/signals.js';
import {divV, span, timestamp} from '../../core/elements.js';
import {text} from '../../core/text.js';
import {handlerRenderer} from '../entities/entity.js';
import {SYSTEM_COLUMNS} from './backend.js';
const REF_ADDRESS = /^\w+\.\w+$/;

/** An action over one row: `requires` names the capability it needs (permission ⇒ hidden),
 * `when` narrows it per row (state ⇒ absent), `run` takes the row. */
export interface DomainAction {
  name: string;
  icon?: string;
  requires?: Capability;
  enabled?: boolean;
  when?: (row: RowView) => boolean;
  run: (row: RowView) => void;
}

export type RowValidator = (value: unknown, row: RowView) => string | null;

export class ActionRegistry {
  private readonly _actions: DomainAction[] = [];

  /** Returns the unregister function. */
  add(action: DomainAction): () => void {
    this._actions.push(action);
    return () => {
      const at = this._actions.indexOf(action);
      if (at >= 0)
        this._actions.splice(at, 1);
    };
  }

  /** The actions that apply to `row`, bound to it — what `rowActions` and a menu take. */
  for(row: RowView): Action[] {
    return ActionRegistry.bind(this._actions, row);
  }

  static bind(actions: readonly DomainAction[], row: RowView): Action[] {
    return actions.filter((a) => a.when === undefined || a.when(row)).map((a): Action =>
      ({name: a.name, icon: a.icon, requires: a.requires, enabled: a.enabled, run: () => a.run(row)}));
  }
}

export class ValidatorRegistry {
  private readonly _byColumn = new Map<string, RowValidator[]>();

  /** Returns the unregister function. */
  add(column: string, validator: RowValidator): () => void {
    const list = this._byColumn.get(column) ?? [];
    this._byColumn.set(column, list);
    list.push(validator);
    return () => {
      const at = list.indexOf(validator);
      if (at >= 0)
        list.splice(at, 1);
    };
  }

  /** The first problem the column's validators report, null when the value passes. */
  check(column: string, value: unknown, row: RowView): string | null {
    for (const validator of this._byColumn.get(column) ?? []) {
      const message = validator(value, row);
      if (message !== null)
        return message;
    }
    return null;
  }
}

export class DomainTable {
  /** Everything an app declares about the table's actions, once. */
  readonly actions = new ActionRegistry();
  readonly validators = new ValidatorRegistry();
  /** The platform handler of the table's rows — deep links, the entity view, `rowFrom`. */
  readonly handler: DG.DomainObjectHandler;
  /** How lists, pickers and chips show a row; the handler's own rendering by default. */
  renderer: ObjectRenderer<RowView>;

  private static readonly _bySource = new WeakMap<DomainSource, DomainTable>();
  private static readonly _byRow = new WeakMap<ReadonlySignal<RowView | null>, DomainSource>();

  constructor(readonly address: string, readonly table: DomainTableLike, readonly access: Access) {
    this.handler = new DG.DomainObjectHandler(address);
    this.renderer = DomainTable.handlerRenderer(this);
  }

  get info(): DomainTableInfoLike {
    return this.table.info;
  }

  get properties(): IProperty[] {
    return this.table.properties;
  }

  /** The handle a source was made through, for the controls over it. */
  static of(source: DomainSource): DomainTable | undefined {
    return DomainTable._bySource.get(source);
  }

  /** The source behind a `currentRow` signal handed to a form. */
  static sourceOf(row: ReadonlySignal<RowView | null>): DomainSource | undefined {
    return DomainTable._byRow.get(row);
  }

  /** Whether the column holds another row's id rather than a value of its own — a ref column
   * (semType `<schema>.<table>`), a user or a group (the js-api editor's `isReferenceProperty`). */
  static isReference(prop: IProperty): boolean {
    const semType = prop.semType ?? '';
    return semType === 'User' || semType === 'Group' || REF_ADDRESS.test(semType);
  }

  /** A started source over this table. */
  source(options: Omit<DomainSourceOptions, 'table'> = {}): DomainSource {
    return this._track(new DomainSource({...options, table: this.address}));
  }

  /** A source holding one pristine draft over `values` as its current row — what a create form
   * binds to; `save()` on it inserts the row. */
  draft(values: Record<string, unknown> = {}): DomainSource {
    return this._track(new DomainSource({table: this.address, draft: true, defaults: values}));
  }

  /** The row as the platform sees it — the handler's `DomainRow`, built locally. */
  row(row: RowView): DG.DomainRow {
    const values: Record<string, unknown> = {...row};
    if (Rows.isDraft(row))
      delete values.id;
    return this.handler.rowFrom(values);
  }

  /** Opens the row's entity view — the platform's default action on a saved row. */
  open(row: RowView): void {
    if (!Rows.isDraft(row))
      this.handler.openRow(this.row(row));
  }

  /** The handler-backed renderer: the name column as the caption where the table declares one,
   * the handler's own rendering otherwise, and the schema-driven rendering where no handler
   * claims the rows. */
  static handlerRenderer(table: DomainTable): ObjectRenderer<RowView> {
    const inner = handlerRenderer<DG.DomainRow>();
    const plain = DomainTable.schemaRenderer(() => ({properties: table.properties, info: table.info}));
    // a draft has no identity the handler could show; a handler card counts only when the handler
    // defines one — the platform's default is the generic property table
    const handled = (row: RowView): DG.DomainRow | null => {
      if (Rows.isDraft(row))
        return null;
      const x = table.row(row);
      return inner.handlerFor(x) === null ? null : x;
    };
    const by = <T>(row: RowView, on: (x: DG.DomainRow) => T, fallback: (row: RowView) => T): T => {
      const x = handled(row);
      return x === null ? fallback(row) : on(x);
    };
    // a card of the handler's own: a registered subclass that overrides renderCard — the reflective
    // default and the platform's Dart meta both paint the generic property table
    const ownCard = (x: DG.DomainRow): boolean => {
      const handler = inner.handlerFor(x);
      if (!(handler instanceof DG.DomainObjectHandler))
        return false;
      for (let p = Object.getPrototypeOf(handler); p !== null && p !== DG.DomainObjectHandler.prototype;
        p = Object.getPrototypeOf(p)) {
        if (Object.prototype.hasOwnProperty.call(p, 'renderCard'))
          return true;
      }
      return false;
    };
    return {
      caption: (row) => {
        const name = plain.caption(row);
        return name === row.id ? by(row, (x) => inner.caption(x), plain.caption) : name;
      },
      icon: (row) => by(row, (x) => inner.icon(x), (r) => span(plain.caption(r))),
      listItem: (row) => by(row, (x) => inner.listItem(x), plain.listItem!),
      markup: (row) => by(row, (x) => inner.markup(x), plain.listItem!),
      tooltip: (row) => by(row, (x) => inner.tooltip(x), plain.card!),
      card: (row) => {
        const x = handled(row);
        return x !== null && ownCard(x) ? inner.card(x) : plain.card!(row);
      },
    };
  }

  /** What a row looks like from its schema alone: the name column (else the first business-key
   * column, else the id; "New <singular>" for a draft without one), and on a card the
   * list-item-rendering recipe — the title, one muted description line, the creation time. */
  static schemaRenderer(schema: () => {properties: IProperty[], info: DomainTableInfoLike}): ObjectRenderer<RowView> {
    const label = (row: RowView): {text: string, draft: boolean} => {
      const info = schema().info;
      const column = info.nameColumn ?? info.businessKey[0];
      const name = column === undefined ? undefined : row[column];
      if (name !== null && name !== undefined && name !== '')
        return {text: String(name), draft: false};
      return Rows.isDraft(row) ? {text: `New ${info.singularName.toLowerCase() || 'row'}`, draft: true} :
        {text: row.id, draft: false};
    };
    const titleOf = (row: RowView): HTMLElement => {
      const {text: caption, draft} = label(row);
      return span(caption, draft ? 'u2-domain-list-name u2-domain-draft' : 'u2-domain-list-name');
    };
    const description = (row: RowView): string => {
      const {properties, info} = schema();
      const prose = properties.filter((p) => p.name !== info.nameColumn &&
        !SYSTEM_COLUMNS.some(([name]) => name === p.name) && !DomainTable.isReference(p) &&
        (p.propertyType ?? p.type) === 'string' && text(row[p.name!]) !== '');
      const first = prose.find((p) => p.name === 'description') ?? prose[0];
      return first === undefined ? '' : text(row[first.name!]);
    };
    return {
      caption: (row) => label(row).text,
      listItem: titleOf,
      card: (row) => {
        const created = row.created_on;
        const line = description(row);
        const title = titleOf(row);
        title.classList.add('u2-domain-card-title');
        return divV([
          title,
          ...(line === '' ? [] : [span(line, 'u2-domain-card-description')]),
          ...(created === null || created === undefined ? [] :
            [timestamp(created as Date | number | string, 'u2-domain-card-time')]),
        ], 'u2-domain-card');
      },
    };
  }

  private _track(source: DomainSource): DomainSource {
    DomainTable._bySource.set(source, this);
    DomainTable._byRow.set(source.currentRow, source);
    source.start();
    return source;
  }
}

export const domains = {
  /** The table's shape and the caller's access in one round of requests; a handle per call, so
   * the access is re-acquired by acquiring again. */
  async table(address: string): Promise<DomainTable> {
    const backend = backends.domain;
    if (backend === undefined)
      throw new Error('no platform backend for domain tables');
    const table = await backend.table(address);
    return new DomainTable(address, table, Access.from(await table.access()));
  },
};
