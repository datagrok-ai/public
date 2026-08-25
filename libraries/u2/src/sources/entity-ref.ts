/* One server entity by id (DD14): the object a form shows the details of. Read-only, so `live` is
   its design-time default; a failed or forbidden lookup ends in the error state and everything
   bound to it reads undefined — the containment posture, never a throw. */
import {signal, computed, Signal} from '../core/signals.js';
import {Component} from '../core/component.js';
import {AsyncValue} from '../core/async-value.js';
import {backends, requireBackend} from './backends.js';
import {AsyncSteps, DesignData} from './async-steps.js';
import {text} from '../core/text.js';
import {subBind} from './sub-bind.js';
import type {BindProp, BindSource} from '../spec/bind-source.js';
import type {ComponentEnv, ComponentStart} from '../spec/registry.js';

export interface EntityRefOptions {
  entityType?: string;
  id?: string;
  designData?: string;
  sample?: unknown;
}

export class EntityRef extends Component implements BindSource, ComponentStart {
  readonly entityType: string;
  /** A signal, so the panel edits it live and a bound path drives it. */
  readonly id: Signal<string>;
  readonly designData: DesignData;
  readonly sample: unknown;

  private readonly _env: ComponentEnv;
  private readonly _async: AsyncValue<unknown>;
  private readonly _status: AsyncSteps<unknown>;
  private readonly _entity: Signal<unknown>;
  private readonly _leaves = new Map<string, Signal<unknown>>();

  constructor(options: EntityRefOptions, env: ComponentEnv) {
    super();
    this._env = env;
    this.entityType = options.entityType ?? 'users';
    this.id = signal(options.id ?? '');
    this.designData = (options.designData as DesignData | undefined) ?? 'live';
    this.sample = options.sample;
    const find = requireBackend(this, backends.dapiFind, 'server entities');

    const live = !env.designTime || this.designData === 'live';
    this._async = new AsyncValue<unknown>(() => this._find(find),
      {auto: live, deps: () => {
        this.id.value;
      }});
    this.own(() => this._async.dispose());
    this._status = new AsyncSteps(this._async, live, this.designData,
      () => this.sample ?? this.componentMeta?.designPreview ?? undefined);
    this._entity = computed(() => {
      const at = this._status.state.value;
      return at.kind === 'ready' ? at.value : undefined;
    }) as unknown as Signal<unknown>;
    this.registerFunction({name: 'refresh', description: 'Look the entity up again', inputs: [],
      apply: () => this._async.refresh()});
  }

  /** Phase two: the id may be bound to an input the form declares after this source. */
  start(): void {
    const bound = subBind(this._env, 'id');
    if (bound !== null)
      this.effect(() => this.id.value = text(bound.value));
  }

  bindStep(name: string): Signal<unknown> | null {
    const status = this._status.step(name);
    if (status !== null)
      return status;
    if (name === 'entity' || name === '')
      return this._entity;
    let leaf = this._leaves.get(name);
    if (leaf === undefined) {
      leaf = computed(() => (this._entity.value as Record<string, unknown> | undefined)?.[name]) as
        unknown as Signal<unknown>;
      this._leaves.set(name, leaf);
    }
    return leaf;
  }

  bindProps(): BindProp[] {
    const props = AsyncSteps.props('Why the last lookup failed');
    props.push({name: 'entity', type: 'object', description: 'The entity itself', default: true});
    const value = this._entity.peek();
    if (value !== null && typeof value === 'object') {
      for (const name of EntityRef._readable(value as object))
        props.push({name, type: null});
    }
    return props;
  }

  private async _find(find: NonNullable<typeof backends.dapiFind>): Promise<unknown> {
    const id = this.id.peek();
    if (id === '')
      throw new Error(`${this.entityType}: no id given`);
    const found = await find(this.entityType, id);
    if (found === null || found === undefined)
      throw new Error(`${this.entityType} "${id}" was not found`);
    return found;
  }

  /** Own keys plus the accessors the value's classes declare: a platform entity keeps every field
   * as a non-enumerable prototype getter, so plain key enumeration would answer almost nothing —
   * and its one own field, `dart`, is the raw handle, not a leaf. */
  private static _readable(value: object): string[] {
    const names = new Set(Object.keys(value));
    for (let at = Object.getPrototypeOf(value); at !== null && at !== Object.prototype;
      at = Object.getPrototypeOf(at)) {
      for (const [name, descriptor] of Object.entries(Object.getOwnPropertyDescriptors(at))) {
        if (descriptor.get !== undefined)
          names.add(name);
      }
    }
    return [...names].filter((name) => name !== 'dart' && !name.startsWith('_'));
  }
}
