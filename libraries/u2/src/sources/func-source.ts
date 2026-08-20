/* A platform function, query or script as data (DD14): the params it is called with are the
   source's binds, its outputs are what the form binds to, and a bound param changing re-runs it.
   Platform-free — everything it needs is the `backends.funcRunner` seam. */
import {signal, computed, Signal, ReadonlySignal} from '../core/signals.js';
import {Component} from '../core/component.js';
import {AsyncValue} from '../core/async-value.js';
import {backends, requireBackend, FuncDescriptorLike} from './backends.js';
import {dfBindings, DataFrameLike} from './df-bindings.js';
import {AsyncSteps, DesignData} from './async-steps.js';
import type {PropertyLike} from '../core/property-like.js';
import type {BindProp, BindSource} from '../spec/bind-source.js';
import type {ComponentEnv, ComponentStart} from '../spec/registry.js';

export type {DesignData} from './async-steps.js';

/** Design-time runs are capped at this many rows — a form is being laid out, not analysed. */
export const DESIGN_ROWS = 100;

const PARAMS = 'params.';

export interface FuncSourceOptions {
  func?: string;
  params?: Record<string, unknown>;
  auto?: boolean;
  debounce?: number;
  designData?: string;
  sample?: unknown;
}

export class FuncSource extends Component implements BindSource, ComponentStart {
  readonly func: string;
  readonly params: Record<string, unknown>;
  readonly auto: boolean;
  readonly debounce: number;
  readonly designData: DesignData;
  readonly sample: unknown;

  private readonly _env: ComponentEnv;
  private readonly _descriptor: FuncDescriptorLike | null;
  private readonly _async: AsyncValue<Record<string, unknown>>;
  private readonly _status: AsyncSteps<Record<string, unknown>>;
  private readonly _bound = signal(new Map<string, Signal<unknown>>());
  private readonly _outputs = new Map<string, Signal<unknown>>();
  private readonly _frames = new Map<string, BindSource>();
  private _sampleValue: Record<string, unknown> | undefined;

  constructor(options: FuncSourceOptions, env: ComponentEnv) {
    super();
    this._env = env;
    this.func = options.func ?? '';
    this.params = {...options.params};
    this.auto = options.auto ?? true;
    this.debounce = options.debounce ?? 300;
    this.sample = options.sample;
    const runner = requireBackend(this, backends.funcRunner, 'functions');
    this._descriptor = this.func === '' ? null : runner.find(this.func);
    const kind = this._descriptor?.kind;
    this.designData = (options.designData as DesignData | undefined) ??
      (kind === 'table-query' || kind === 'query' ? 'live' : 'schema');

    const live = !env.designTime || this.designData === 'live';
    // DD9's letter: at design time a function that is not known to READ is never AUTO-run, however
    // live the policy — the tray's explicit Refresh is the only run it gets. Only a visual query
    // is known to read: a `DataQuery` carries hand-written SQL the connector executes as written,
    // DML included, and opening someone's spec must not run it
    const kicks = this.auto && live && (!env.designTime || kind === 'table-query');
    this._async = new AsyncValue<Record<string, unknown>>((abort) => this._run(runner, abort),
      {auto: kicks, debounceMs: this.debounce, deps: () => this._deps()});
    this.own(() => this._async.dispose());
    this._status = new AsyncSteps(this._async, live, this.designData, () => this._sample());
    this.registerFunction({name: 'refresh', description: 'Run the function again', inputs: [],
      apply: () => this._async.refresh()});
  }

  /** Phase two: the form exists, so a param bound to one of its inputs resolves now. */
  start(): void {
    const bound = new Map<string, Signal<unknown>>();
    for (const [key, path] of Object.entries(this._env.subBinds)) {
      if (!key.startsWith(PARAMS))
        continue;
      const source = this._env.resolve(path);
      if (source !== null)
        bound.set(key.slice(PARAMS.length), source);
    }
    if (bound.size > 0)
      this._bound.value = bound;
  }

  bindStep(name: string): Signal<unknown> | BindSource | null {
    const status = this._status.step(name);
    if (status !== null)
      return status;
    const key = name === '' ? this._defaultOutput() : name;
    if (key !== null && this._declares(key))
      return this._isFrame(key) ? this._frame(key) : this._output(key);
    // a single-table source IS that table: `$.orders.currentRow.name` skips the output's name,
    // which the author of the form neither chose nor cares about
    const single = this._defaultOutput();
    return name !== '' && single !== null && this._isFrame(single) ?
      this._frame(single).bindStep(name) : null;
  }

  /** Read field by field, never spread: a platform `DG.Property` keeps everything on getters, so a
   * copy of it is an empty object with a `dart` handle — and an output with no name is a blank row
   * in the binding picker. */
  bindProps(): BindProp[] {
    const props = AsyncSteps.props('Why the last run failed');
    for (const output of this._descriptor?.outputs ?? []) {
      props.push({name: output.name, type: output.propertyType ?? output.type,
        semType: output.semType, description: output.description,
        walkable: FuncSource._isFrameProp(output)});
    }
    return props;
  }

  private _run(runner: NonNullable<typeof backends.funcRunner>,
    abort: AbortSignal): Promise<Record<string, unknown>> {
    if (this._descriptor === null)
      return Promise.reject(new Error(`"${this.func}" is not a known function`));
    const values = {...this.params};
    for (const [name, source] of this._bound.peek())
      values[name] = source.peek();
    return runner.run(this.func, values, abort,
      this._env.designTime ? {maxRows: DESIGN_ROWS} : undefined);
  }

  /** Read inside the AsyncValue's effect: every bound param is a dependency, so an upstream edit
   * cancels the run in flight and starts a new one. */
  private _deps(): void {
    for (const source of this._bound.value.values())
      source.value;
  }

  private _defaultOutput(): string | null {
    const outputs = this._descriptor?.outputs ?? [];
    return outputs.length === 1 ? outputs[0].name : null;
  }

  /** A known function offers exactly what it declares; an unknown one takes any name and answers
   * undefined, so a bind survives until the function is chosen. */
  private _declares(name: string): boolean {
    return this._descriptor === null ||
      this._descriptor.outputs.some((output) => output.name === name);
  }

  private _isFrame(name: string): boolean {
    const output = this._descriptor?.outputs.find((p) => p.name === name);
    return output !== undefined && FuncSource._isFrameProp(output);
  }

  private _output(name: string): Signal<unknown> {
    let output = this._outputs.get(name);
    if (output === undefined) {
      const state = this._status.state;
      output = computed(() => {
        const at = state.value;
        return at.kind === 'ready' ? at.value[name] : undefined;
      }) as unknown as Signal<unknown>;
      this._outputs.set(name, output);
    }
    return output;
  }

  private _frame(name: string): BindSource {
    let frame = this._frames.get(name);
    if (frame === undefined) {
      const source = this._output(name) as unknown as ReadonlySignal<DataFrameLike | undefined>;
      frame = dfBindings(source, this.scope);
      this._frames.set(name, frame);
    }
    return frame;
  }

  /** Built once and kept: a fresh frame on every read would repoint every binding under it, and
   * the fields the preview exists to fill would blink back to empty. */
  private _sample(): Record<string, unknown> {
    if (this._sampleValue === undefined)
      this._sampleValue = this._sampleOutputs();
    return this._sampleValue;
  }

  /** The sample stored with the spec (C5), or the tag's own fallback: output name → value, or a
   * bare array of rows for a single-output function. Rows become a frame, so the canvas previews
   * exactly what a live run would show. */
  private _sampleOutputs(): Record<string, unknown> {
    const sample = this.sample ?? this.meta?.designPreview;
    if (sample === null || sample === undefined)
      return {};
    const single = this._defaultOutput();
    const record = Array.isArray(sample) ?
      (single === null ? {} : {[single]: sample}) : sample as Record<string, unknown>;
    const outputs: Record<string, unknown> = {};
    for (const [name, value] of Object.entries(record))
      outputs[name] = FuncSource._rows(value) ?? value;
    return outputs;
  }

  private static _rows(value: unknown): DataFrameLike | undefined {
    const make = backends.tableFromRows;
    return make !== undefined && Array.isArray(value) && value.length > 0 &&
      value.every((row) => typeof row === 'object' && row !== null) ?
      make(value as object[]) : undefined;
  }

  private static _isFrameProp(prop: PropertyLike): boolean {
    return (prop.propertyType ?? prop.type) === 'dataframe';
  }
}
