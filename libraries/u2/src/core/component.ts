import {Scope} from './scope.js';
import {Signal, computed, isWritableSignal, signal} from './signals.js';
import {propertyFields} from './property-like.js';
import type {PropertyLike} from './property-like.js';
import type {BindPropLike, BindSourceLike, ComponentMetaLike, FuncLike, ObservableLike, WidgetLike,
  WidgetStatusLike} from './widget-like.js';

/** What a widget announces when one of its properties changed: the property's name, or null when
 * several changed at once (or the widget cannot say which). Structural — nothing platform-typed. */
export type PropertyChange = string | null;

/** One property-tier step: `inner` is what the step reads, `known` the value the tier last
 * observed — read in a refresh or written by the effect — so the write effect never compares
 * against a fresh read (a read answering a new wrapper per call would spin otherwise). */
interface PropertyEntry {
  inner: Signal<unknown>;
  step: Signal<unknown>;
  known: unknown;
}

/** The per-instance u2 state, kept under one own field so that {@link Component.init} can be
 * applied to a foreign object (an adopted `DG.Widget`) without clobbering its own fields. */
export interface ComponentState {
  meta?: ComponentMetaLike;
  aiDescription: string | null;
  properties: PropertyLike[] | null;
  propertiesMeta?: ComponentMetaLike;
  functions: FuncLike[];
  listeners: {id: string | null, next: (x: unknown) => void}[];
  propertySignals: Map<string, PropertyEntry>;
  propertyIndex?: Map<string, PropertyLike>;
  warned?: Set<string>;
  changesWired: boolean;
  /** Names a notification marked since the last refresh; 'all' for an anonymous change. */
  stale?: Set<string> | 'all';
  table?: unknown;
  detaching?: boolean;
  platformKilled?: boolean;
}

/** Base of every u2 component, visual or not: a name, an effect/cleanup scope, and the widget
 * introspection surface (DD7) generated from {@link meta}. A component constructed inside
 * another's scope is disposed with it. Platform-free; in Datagrok, `u2/dg`'s `host()` wires
 * {@link dispose} into the existing `DG.Widget` kill channel and delegates the introspection to a
 * real `DG.Widget`. Standalone hosts (gallery, tests) call {@link dispose} directly. */
export class Component implements WidgetLike, BindSourceLike {
  readonly scope!: Scope;
  /** Unique within a spec: the anchor for selection, patches and automation. */
  name?: string;
  /** The registry metadata of the tag that built this component, stamped by the spec renderer —
   * the single source the introspection surface below is generated from. */
  meta?: ComponentMetaLike;
  /** The props the registry's `create` was handed, stamped alongside {@link meta}: the last resort
   * of the property read below, for a component that keeps a prop under another name or not at all
   * (a `create` wrapping a bare element). Construction-time values — a live one always wins. */
  specProps?: Record<string, unknown>;
  readonly _u2!: ComponentState;

  constructor(scope?: Scope) {
    Component.init(this, scope);
  }

  /** Per-instance state, for a constructed component and for an adopted platform object alike. */
  static init(self: Component, scope?: Scope): void {
    const it = self as {scope: Scope, _u2: ComponentState};
    it.scope = scope ?? new Scope();
    it._u2 = {aiDescription: null, properties: null, functions: [], listeners: [],
      propertySignals: new Map(), changesWired: false};
    const owner = Scope.ambient;
    if (owner && owner !== it.scope)
      owner.own(() => self.dispose());
  }

  static is(x: unknown): x is Component {
    if (x instanceof Component)
      return true;
    const c = x as Partial<Component> | null | undefined;
    return typeof c?.bindStep === 'function' && typeof c?.bindProps === 'function' && c?.scope instanceof Scope;
  }

  /** The equality the property tier and {@link link} suppress echoes by: `Object.is`, arrays
   * element-wise, plain objects key-wise (recursive), everything else by identity. */
  static sameValue(a: unknown, b: unknown): boolean {
    if (Object.is(a, b))
      return true;
    if (Array.isArray(a) && Array.isArray(b))
      return a.length === b.length && a.every((x, i) => Component.sameValue(x, b[i]));
    if (!Component._isPlain(a) || !Component._isPlain(b))
      return false;
    const keys = Object.keys(a);
    return keys.length === Object.keys(b).length &&
      keys.every((k) => k in b && Component.sameValue(a[k], b[k]));
  }

  private static _isPlain(x: unknown): x is Record<string, unknown> {
    if (x === null || typeof x !== 'object')
      return false;
    const proto = Object.getPrototypeOf(x);
    return proto === Object.prototype || proto === null;
  }

  run<T>(fn: () => T): T {
    return Scope.runWith(this.scope, fn);
  }

  effect(fn: () => void): this {
    this.scope.effect(fn);
    return this;
  }

  own(dispose: () => void): this {
    this.scope.own(dispose);
    return this;
  }

  dispose(): void {
    this.scope.dispose();
  }

  /** A short AI-facing briefing, seeded from the registry's usage or description. */
  get aiDescription(): string | null {
    return this._u2.aiDescription ?? this.meta?.usage ?? this.meta?.description ?? null;
  }

  set aiDescription(x: string | null) {
    this._u2.aiDescription = x;
  }

  /** The props {@link meta} declares, each closed over the component's own state: a same-named
   * signal reads and writes live, a same-named plain member, an `options` entry or — last —
   * a {@link specProps} entry reads read-only. {@link name} — the component's spec identity —
   * is never a property; an input's
   * declared `name` prop is its form key, a different thing that IS a property like any other.
   * A component never hand-writes its property list. The list is generated once per {@link meta}:
   * the accessors read live state, so only a re-stamp invalidates it. */
  getProperties(): PropertyLike[] {
    const u2 = this._u2;
    if (u2.properties && u2.propertiesMeta === this.meta)
      return u2.properties;
    const self = this as unknown as Record<string, unknown>;
    u2.propertiesMeta = this.meta;
    u2.propertyIndex = undefined;
    u2.properties = (this.meta?.props ?? []).map((p) => {
      const prop: PropertyLike = {...p, get: () => Component._read(self, p.name)};
      if (self[p.name] instanceof Signal) {
        prop.set = (_source: unknown, value: unknown) => {
          (self[p.name] as Signal<unknown>).value = value;
        };
      }
      return prop;
    });
    return u2.properties;
  }

  getFunctions(): FuncLike[] {
    return this._u2.functions;
  }

  /** Events fired by this component, and — on a {@link Control} — the DOM events of that name
   * raised on its root. Without an id, every component event. The two streams are merged as they
   * are, so a subscriber sees either payload and must discriminate (see {@link ObservableLike}). */
  onEvent(eventId: string | null = null): ObservableLike<unknown> {
    return {
      subscribe: (next: (x: unknown) => void) => {
        const listeners = this._u2.listeners;
        const listener = {id: eventId, next};
        listeners.push(listener);
        const off = () => {
          const i = listeners.indexOf(listener);
          if (i >= 0)
            listeners.splice(i, 1);
        };
        this.scope.own(off);
        return {unsubscribe: () => {
          off();
          this.scope.disown(off);
        }};
      },
    };
  }

  getWidgetStatus(): WidgetStatusLike {
    return {
      parts: {},
      hitAreas: {},
      shortcuts: {},
      events: (this.meta?.events ?? []).map((name) => ({name, eventName: name, description: ''})),
      description: this.meta?.description ?? null,
      error: null,
    };
  }

  /** The second binding tier is on: every property {@link getProperties} declares is a bind step,
   * read and written on {@link propertyTarget}. Off for a plain u2 component, whose signal-less
   * meta props are construction-time options. */
  get propertyTier(): boolean {
    return false;
  }

  /** What this component's properties are read from and written to — `prop.get(target)` /
   * `prop.set(target, v)`: the component itself, unless an adoption says otherwise. */
  get propertyTarget(): unknown {
    return this;
  }

  /** What refreshes the property tier's cached signals from outside: null where nothing announces
   * a change (one-way). */
  protected propertyChanges(): ObservableLike<PropertyChange> | null {
    return null;
  }

  /** Keeps this component's `name` step and `source` equal, echo-suppressed by value; the reverse
   * direction only where `twoWay`. A component that took `source` as its own signal needs no link
   * and gets none; a read-only step cannot follow anything and gets none either. */
  link(name: string, source: Signal<unknown>, twoWay: boolean): void {
    const own = this.bindStep(name);
    if (!(own instanceof Signal) || own === source)
      return;
    if (!isWritableSignal(own)) {
      const owner = this.name ?? this.meta?.tag ?? 'component';
      this._warnOnce(name, `${owner}.${name} is read-only — edits will not flow back`);
      return;
    }
    const scope = Scope.ambient ?? this.scope;
    scope.effect(() => {
      const value = source.value;
      if (!Component.sameValue(own.peek(), value))
        own.value = value;
    });
    if (!twoWay)
      return;
    scope.effect(() => {
      const value = own.value;
      if (!Component.sameValue(source.peek(), value))
        source.value = value;
    });
  }

  /** One binding step (DD5): a meta-declared prop's same-named Signal member — a named node is a
   * bind source, generated off {@link meta} the way {@link getProperties} is — then, under
   * {@link propertyTier}, a cached signal over the declared property. '' answers the default
   * binding: the component's own `value` signal, where one exists. */
  bindStep(name: string): Signal<unknown> | BindSourceLike | null {
    const self = this as unknown as Record<string, unknown>;
    const member = name === '' ? self.value :
      this.meta?.props?.some((p) => p.name === name) ? self[name] : undefined;
    if (member instanceof Signal)
      return member as Signal<unknown>;
    return this.propertyTier && name !== '' ? this._propertySignal(name) : null;
  }

  /** The signal-backed subset of the meta props — what a binding picker offers on a named node —
   * plus, under {@link propertyTier}, every declared property. Allocates no signal. */
  bindProps(): BindPropLike[] {
    const self = this as unknown as Record<string, unknown>;
    const props: BindPropLike[] = (this.meta?.props ?? [])
      .filter((p) => self[p.name] instanceof Signal)
      .map((p) => ({...p, writable: isWritableSignal(self[p.name]),
        ...(p.name === 'value' ? {default: true} : {})}));
    if (!this.propertyTier)
      return props;
    const own = new Set(props.map((p) => p.name));
    for (const p of this._propertyIndex().values()) {
      if (!own.has(p.name))
        props.push({...propertyFields(p), writable: p.set != null});
    }
    return props;
  }

  private _property(name: string): PropertyLike | undefined {
    return this.propertyTier ? this._propertyIndex().get(name) :
      this.getProperties().find((p) => p.name === name);
  }

  private _propertyIndex(): Map<string, PropertyLike> {
    const u2 = this._u2;
    if (u2.propertyIndex === undefined)
      u2.propertyIndex = new Map(this.getProperties().map((p) => [p.name, p]));
    return u2.propertyIndex;
  }

  private _propertySignal(name: string): Signal<unknown> | null {
    const u2 = this._u2;
    const cached = u2.propertySignals.get(name);
    if (cached)
      return cached.step;
    const prop = this._property(name);
    if (prop === undefined)
      return null;
    const entry = {known: prop.get?.(this.propertyTarget)} as PropertyEntry;
    entry.inner = signal(entry.known);
    if (prop.set == null)
      entry.step = computed(() => entry.inner.value) as unknown as Signal<unknown>;
    else {
      entry.step = entry.inner;
      this.scope.effect(() => {
        const value = entry.inner.value;
        if (Component.sameValue(entry.known, value))
          return;
        entry.known = value;
        prop.set!(this.propertyTarget, value);
      });
    }
    u2.propertySignals.set(name, entry);
    this._wireChanges();
    return entry.step;
  }

  private _wireChanges(): void {
    const u2 = this._u2;
    if (u2.changesWired)
      return;
    u2.changesWired = true;
    const changes = this.propertyChanges();
    if (changes === null)
      return;
    const subscription = changes.subscribe((name) => this._stale(name));
    this.scope.own(() => subscription.unsubscribe());
  }

  /** A notification never refreshes inside the event: the names pile up and one microtask reads
   * them back — the platform fires on every write, the effect's own included. */
  private _stale(name: PropertyChange): void {
    const u2 = this._u2;
    const first = u2.stale === undefined;
    if (name === null)
      u2.stale = 'all';
    else if (u2.stale === undefined)
      u2.stale = new Set([name]);
    else if (u2.stale !== 'all')
      u2.stale.add(name);
    if (first)
      queueMicrotask(() => this._refreshStale());
  }

  private _refreshStale(): void {
    const u2 = this._u2;
    const stale = u2.stale;
    u2.stale = undefined;
    if (stale === undefined || this.scope.isDisposed)
      return;
    for (const name of stale === 'all' ? u2.propertySignals.keys() : stale)
      this._refreshProperty(name);
  }

  private _refreshProperty(name: string): void {
    const entry = this._u2.propertySignals.get(name);
    if (entry === undefined)
      return;
    const value = this._property(name)?.get?.(this.propertyTarget);
    entry.known = value;
    if (!Component.sameValue(entry.inner.peek(), value))
      entry.inner.value = value;
  }

  private _warnOnce(key: string, message: string): void {
    const warned = this._u2.warned ??= new Set();
    if (warned.has(key))
      return;
    warned.add(key);
    console.warn(`u2: ${message}`);
  }

  private static _read(self: Record<string, unknown>, name: string): unknown {
    const member = self[name];
    if (member instanceof Signal)
      return member.peek();
    if (member !== undefined)
      return member;
    const option = (self.options as Record<string, unknown> | undefined)?.[name];
    if (option !== undefined)
      return option;
    const given = (self.specProps as Record<string, unknown> | undefined)?.[name];
    return given instanceof Signal ? given.peek() : given;
  }

  protected registerFunction(f: FuncLike): void {
    this._u2.functions.push(f);
  }

  protected fireEvent(eventId: string, args?: unknown): void {
    for (const listener of this._u2.listeners.slice()) {
      if (listener.id === null || listener.id === eventId)
        listener.next(args);
    }
  }
}

/** A component that renders: everything visual in u2 is one. */
export class Control extends Component {
  readonly root: HTMLElement;

  constructor(root?: HTMLElement, scope?: Scope) {
    super(scope);
    this.root = root ?? document.createElement('div');
    this.root.classList.add('u2-component');
  }

  static is(x: unknown): x is Control {
    return (x as Partial<Control> | null | undefined)?.root instanceof HTMLElement && Component.is(x);
  }

  /** Component events plus the DOM events of that name raised on {@link root} — the listener is
   * wired per subscription and dropped with it, or with the scope if the subscriber never does. */
  onEvent(eventId: string | null = null): ObservableLike<unknown> {
    const events = super.onEvent(eventId);
    if (eventId === null)
      return events;
    return {
      subscribe: (next: (x: unknown) => void) => {
        const subscription = events.subscribe(next);
        const handler = (e: Event) => next(e);
        this.root.addEventListener(eventId, handler);
        const off = () => this.root.removeEventListener(eventId, handler);
        this.scope.own(off);
        return {unsubscribe: () => {
          subscription.unsubscribe();
          off();
          this.scope.disown(off);
        }};
      },
    };
  }

  /** Builds a control out of plain elements: `builder` runs with the new control's scope
   * ambient, so signal bindings and nested components inside it are owned by the result.
   * The ambient scope only spans the synchronous part of `builder` — components constructed after
   * an `await` inside it are adopted by nobody and must be disposed by hand, or constructed before
   * the first await. */
  static build(builder: () => HTMLElement | (HTMLElement | Control | string)[]): Control {
    const scope = new Scope();
    const built = Scope.runWith(scope, builder);
    if (!Array.isArray(built))
      return new Control(built, scope);
    const root = document.createElement('div');
    for (const child of built)
      root.append(Control.is(child) ? child.root : child);
    return new Control(root, scope);
  }
}
