/* An open workspace table as data (DD14): the whole dfBindings surface over the frame the
   platform currently answers to that name, following the workspace as tables open and close.
   A name is a weak reference by design (DATA_BINDING §8) — an absent table warns and the bound
   fields stay empty rather than breaking. */
import {signal, Signal} from '../core/signals.js';
import {Component} from '../core/component.js';
import {specWarn} from '../spec/spec.js';
import {backends, requireBackend} from './backends.js';
import {dfBindings, DataFrameLike} from './df-bindings.js';
import type {BindProp, BindSource} from '../spec/bind-source.js';

export interface TableSourceOptions {
  table?: string;
}

export class TableSource extends Component implements BindSource {
  readonly table: string;

  private readonly _df: Signal<DataFrameLike | undefined>;
  private readonly _bindings: BindSource;

  constructor(options: TableSourceOptions = {}) {
    super();
    this.table = options.table ?? '';
    const workspace = requireBackend(this, backends.workspace, 'the workspace');
    const read = () => workspace.table(this.table) ?? undefined;
    this._df = signal(read());
    const sub = workspace.onTablesChanged.subscribe(() => this._df.value = read());
    this.own(() => sub.unsubscribe());
    if (this._df.peek() === undefined)
      specWarn(`u2-table-source: no open table named "${this.table}"`);
    this._bindings = dfBindings(this._df, this.scope);
  }

  bindStep(name: string): Signal<unknown> | BindSource | null {
    return this._bindings.bindStep(name);
  }

  bindProps(): BindProp[] {
    return this._bindings.bindProps();
  }
}
