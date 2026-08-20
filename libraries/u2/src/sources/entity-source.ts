/* A server collection as data (DD14): the platform's own paging over `grok.dapi.*`, exposed as
   the items, their labels, the paging state and the total. Always live — a dapi listing is
   read-only and `pageSize` is its own cap, so the design-time policies do not apply (Q5). */
import {signal, computed, Signal} from '../core/signals.js';
import {Component} from '../core/component.js';
import {AsyncPager} from '../core/pager.js';
import {dapiPager, sanitizeFilterValue} from '../dg/entities/dapi-source.js';
import {backends, requireBackend} from './backends.js';
import {text} from '../core/text.js';
import {subBind} from './sub-bind.js';
import type {BindProp, BindSource} from '../spec/bind-source.js';
import type {ComponentEnv, ComponentStart} from '../spec/registry.js';

export interface EntitySourceOptions {
  entity?: string;
  filter?: string;
  order?: string;
  pageSize?: number;
}

export class EntitySource extends Component implements BindSource, ComponentStart {
  readonly entity: string;
  readonly order: string | undefined;
  readonly pageSize: number;
  /** A signal, so the panel edits it live and a bound path drives it. */
  readonly filter: Signal<string>;

  private readonly _env: ComponentEnv;
  private readonly _pager: AsyncPager<unknown>;
  private readonly _steps: Record<string, Signal<unknown>>;

  constructor(options: EntitySourceOptions, env: ComponentEnv) {
    super();
    this._env = env;
    this.entity = options.entity ?? 'users';
    this.order = options.order;
    this.pageSize = options.pageSize ?? 20;
    this.filter = signal(options.filter ?? '');
    const collection = requireBackend(this, backends.dapi, 'server collections');
    this._pager = dapiPager<unknown>(() => {
      const source = collection(this.entity);
      if (source === null)
        throw new Error(`unknown collection "${this.entity}"`);
      return source;
    }, {pageSize: this.pageSize, order: this.order,
      filter: () => this.filter.peek() === '' ? undefined : this.filter.peek()});
    this.own(() => this._pager.dispose());

    const items = this._pager.items as unknown as Signal<unknown>;
    this._steps = {
      items,
      names: computed(() => this._pager.items.value.map(EntitySource._label)) as unknown as Signal<unknown>,
      state: this._pager.state as unknown as Signal<unknown>,
      total: this._pager.total as unknown as Signal<unknown>,
    };
    this.registerFunction({name: 'refresh', description: 'Load the collection again', inputs: [],
      apply: () => this._pager.reset()});
    this.registerFunction({name: 'loadMore', description: 'Append the next page', inputs: [],
      apply: () => this._pager.loadMore()});
  }

  /** Phase two: the filter may be bound to an input the form declares after this source, and the
   * first page is loaded once that value is known. A bound value is whatever the user typed into
   * that input, so it is sanitized before it reaches the smart-filter grammar; the literal `filter`
   * prop stays raw — the author of the spec writes that one, and writes it in the grammar. */
  start(): void {
    const bound = subBind(this._env, 'filter');
    if (bound !== null)
      this.effect(() => this.filter.value = sanitizeFilterValue(text(bound.value)));
    // a new filter is a new collection: reset drops the accumulated pages and reloads from page 0
    this.effect(() => {
      this.filter.value;
      this._pager.reset();
    });
  }

  bindStep(name: string): Signal<unknown> | null {
    return this._steps[name === '' ? 'items' : name] ?? null;
  }

  bindProps(): BindProp[] {
    return [
      {name: 'items', type: 'object', description: 'The entities loaded so far'},
      {name: 'names', type: 'string_list', description: 'Their names, for a choice input'},
      {name: 'state', type: 'string', description: 'idle, loading, done or error'},
      {name: 'total', type: 'int'},
    ];
  }

  private static _label(item: unknown): string {
    const entity = item as {friendlyName?: string, name?: string} | null;
    return entity?.friendlyName ?? entity?.name ?? String(item);
  }
}
