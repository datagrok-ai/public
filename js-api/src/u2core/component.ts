import {Scope} from './scope.js';
import {Signal, computed, isWritableSignal, signal} from './signals.js';
import {propertyFields} from './property-fields.js';
import type {IProperty} from './property-fields.js';
import type {BindProp, BindSource, ComponentMetaBase, FuncLike, NamedProperty, ObservableLike,
  IWidgetStatus} from './protocol.js';

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
 * applied to a platform `DG.Widget` (a Component by inheritance) without clobbering its own fields. */
export interface ComponentState {
  scope?: Scope;
  aiDescription: string | null;
  properties: NamedProperty[] | null;
  propertiesMeta?: ComponentMetaBase;
  functions: FuncLike[];
  listeners: {id: string | null, next: (x: unknown) => void}[];
  propertySignals: Map<string, PropertyEntry>;
  propertyIndex?: Map<string, NamedProperty>;
  warned?: Set<string>;
  changesWired: boolean;
  /** Names a notification marked since the last refresh; 'all' for an anonymous change. */
  stale?: Set<string> | 'all';
  table?: unknown;
  detaching?: boolean;
  platformKilled?: boolean;
  propertyChangeListeners?: Set<(c: PropertyChange) => void>;
  name?: string;
}

/** Base of every u2 component, visual or not: a name, an effect/cleanup scope, and the widget
 * introspection surface (DD7) generated from {@link componentMeta} — what the shell, the context
 * panel, automation and the copilot ask of a widget, implemented by `DG.Widget` through
 * inheritance. A component constructed inside
 * another's scope is disposed with it. Platform-free; in Datagrok, `DG.Widget extends Control`,
 * so the kill channel and the introspection surface are this class's own — nothing wires or
 * delegates. Standalone hosts (gallery, tests) call {@link dispose} directly. */
export class Component implements BindSource {
  /** Minted on first use; assigned before {@link _scopeMinted} fires, so the hook's own
   * `this.scope.own(...)` reads the stored scope instead of re-entering the mint. */
  get scope(): Scope {
    const u2 = this._u2;
    if (u2.scope === undefined) {
      u2.scope = new Scope();
      this._scopeMinted(u2.scope);
    }
    return u2.scope;
  }

  /** First-mint hook: a platform descendant wires its lifecycle here. Never fired for a
   * scope provided to {@link init} — that path runs inside the constructor chain. */
  protected _scopeMinted(_scope: Scope): void {}
  /** Unique within a spec: the anchor for selection, patches and automation. An accessor pair
   * (not a field), so a platform descendant's own `name` accessors override it legally (TS2611). */
  get name(): string | undefined { return this._u2.name; }
  set name(x: string | undefined) { this._u2.name = x; }

  /** Sets the introspection/form key without claiming an automation id — the label-fallback path:
   * a caption is not a stable identity, so it never reaches `data-u2-name`. */
  protected deriveName(x: string | undefined): void {
    this._u2.name = x;
  }
  /** The registry metadata of the tag that built this component, stamped by the spec renderer —
   * the single source the introspection surface below is generated from. */
  componentMeta?: ComponentMetaBase;
  /** The props the registry's `create` was handed, stamped alongside {@link componentMeta}: the last resort
   * of the property read below, for a component that keeps a prop under another name or not at all
   * (a `create` wrapping a bare element). Construction-time values — a live one always wins. */
  specProps?: Record<string, unknown>;
  readonly _u2!: ComponentState;

  constructor(scope?: Scope) {
    Component.init(this, scope);
  }

  /** Per-instance state, for a constructed component and for an adopted platform object alike. */
  static init(self: Component, scope?: Scope): void {
    const it = self as {_u2: ComponentState};
    it._u2 = {aiDescription: null, properties: null, functions: [], listeners: [],
      propertySignals: new Map(), changesWired: false};
    if (scope !== undefined)
      it._u2.scope = scope;
    const owner = self.ambientAdoption ? Scope.ambient : undefined;
    if (owner && owner !== it._u2.scope)
      owner.own(() => self.dispose());
  }

  /** Structural probe by methods only — never `.scope`: the accessor mints. */
  static is(x: unknown): x is Component {
    if (x instanceof Component)
      return true;
    const c = x as Partial<Component> | null | undefined;
    return typeof c?.bindStep === 'function' && typeof c?.bindProps === 'function' &&
      typeof c?.own === 'function' && typeof c?.dispose === 'function';
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

  runInScope<T>(fn: () => T): T {
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

  /** On a never-engaged wrapper the getter mints here and {@link _scopeMinted} wires the
   * lifecycle (`own(kill)` + the platform cleanup for Dart-owned), so dispose still kills it;
   * re-entry from the kill cleanup is a no-op via the scope's disposed flag. */
  dispose(): void {
    this.scope.dispose();
  }

  /** A short AI-facing briefing, seeded from the registry's usage or description. */
  get aiDescription(): string | null {
    return this._u2.aiDescription ?? this.componentMeta?.usage ?? this.componentMeta?.description ?? null;
  }

  set aiDescription(x: string | null) {
    this._u2.aiDescription = x;
  }

  /** The props {@link componentMeta} declares, each closed over the component's own state: a same-named
   * signal reads and writes live, a same-named plain member, an `options` entry or — last —
   * a {@link specProps} entry reads read-only. {@link name} — the component's spec identity —
   * is never a property; an input's
   * declared `name` prop is its form key, a different thing that IS a property like any other.
   * A component never hand-writes its property list. The list is generated once per {@link componentMeta}:
   * the accessors read live state, so only a re-stamp invalidates it. */
  getProperties(): NamedProperty[] {
    const u2 = this._u2;
    if (u2.properties && u2.propertiesMeta === this.componentMeta)
      return u2.properties;
    const self = this as unknown as Record<string, unknown>;
    u2.propertiesMeta = this.componentMeta;
    u2.propertyIndex = undefined;
    u2.properties = (this.componentMeta?.props ?? []).map((p) => {
      const name = p.name;
      const prop: NamedProperty = {...p, get: () => Component._read(self, name)};
      if (self[name] instanceof Signal) {
        prop.set = (_source: unknown, value: unknown) => {
          (self[name] as Signal<unknown>).value = value;
        };
      } else if (Component._settable(self, name)) {
        prop.set = (_source: unknown, value: unknown) => {
          self[name] = value;
        };
      }
      return prop;
    });
    return u2.properties;
  }

  /** A declared prop backed by an accessor with a setter is writable through it — the setter
   * carries the live behavior. A plain data field stays read-only: a bare assignment would change
   * nothing the component watches. */
  private static _settable(self: object, name: string): boolean {
    for (let at: object | null = self; at !== null; at = Object.getPrototypeOf(at)) {
      const descriptor = Object.getOwnPropertyDescriptor(at, name);
      if (descriptor !== undefined)
        return descriptor.set !== undefined;
    }
    return false;
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

  getWidgetStatus(): IWidgetStatus {
    return {
      parts: {},
      hitAreas: {},
      shortcuts: {},
      events: (this.componentMeta?.events ?? []).map((name) => ({name, eventName: name, description: ''})),
      description: this.componentMeta?.description ?? null,
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

  /** Whether a component constructed inside an ambient scope is disposed with it. */
  protected get ambientAdoption(): boolean {
    return true;
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
      const owner = this.name ?? this.componentMeta?.tag ?? 'component';
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
   * bind source, generated off {@link componentMeta} the way {@link getProperties} is — then, under
   * {@link propertyTier}, a cached signal over the declared property. '' answers the default
   * binding: the component's own `value` signal, where one exists. */
  bindStep(name: string): Signal<unknown> | BindSource | null {
    const self = this as unknown as Record<string, unknown>;
    const member = name === '' ? self.value :
      this.componentMeta?.props?.some((p) => p.name === name) ? self[name] : undefined;
    if (member instanceof Signal)
      return member as Signal<unknown>;
    return this.propertyTier && name !== '' ? this._propertySignal(name) : null;
  }

  /** The signal-backed subset of the meta props — what a binding picker offers on a named node —
   * plus, under {@link propertyTier}, every declared property. Allocates no signal. */
  bindProps(): BindProp[] {
    const self = this as unknown as Record<string, unknown>;
    const props: BindProp[] = (this.componentMeta?.props ?? [])
      .filter((p) => self[p.name] instanceof Signal)
      .map((p) => ({...p, writable: isWritableSignal(self[p.name]),
        ...(p.name === 'value' ? {default: true} : {})}));
    if (!this.propertyTier)
      return props;
    const own = new Set(props.map((p) => p.name));
    for (const p of this._propertyIndex().values()) {
      if (!own.has(p.name))
        props.push({...propertyFields(p), name: p.name, writable: p.set != null});
    }
    return props;
  }

  private _property(name: string): NamedProperty | undefined {
    return this.propertyTier ? this._propertyIndex().get(name) :
      this.getProperties().find((p) => p.name === name);
  }

  private _propertyIndex(): Map<string, NamedProperty> {
    const u2 = this._u2;
    if (u2.propertyIndex === undefined)
      u2.propertyIndex = new Map(this.getProperties().map((p): [string, NamedProperty] => [p.name, p]));
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
    const scope = u2.scope;
    if (stale === undefined || scope === undefined || scope.isDisposed)
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
  protected _root!: HTMLElement;
  private static readonly _byRoot = new WeakMap<Element, Control>();

  constructor(root?: HTMLElement, scope?: Scope) {
    super(scope);
    this._root = root ?? document.createElement('div');
    this._root.classList.add('u2-component');
    Control._byRoot.set(this._root, this);
  }

  get root(): HTMLElement {
    return this._root;
  }

  set root(x: HTMLElement) {
    this._root = x;
    Control._byRoot.set(x, this);
    this._stampName();
  }

  /** A named control carries its name in the DOM: `data-u2-name` is the stable automation id
   * selectors, tutorials and agents address the control by — the same attribute the spec renderer
   * stamps on named nodes, so hand-built and spec-built UI share one selector vocabulary.
   * Derived names (an input's label fallback) never stamp — see {@link Component.deriveName}. */
  get name(): string | undefined {
    return this._u2.name;
  }

  set name(x: string | undefined) {
    this._u2.name = x;
    this._stampName();
  }

  private _stampName(): void {
    const name = this._u2.name;
    if (name === undefined)
      delete this._root.dataset.u2Name;
    else
      this._root.dataset.u2Name = name;
  }

  /** The nearest control owning `el` — the DOM → component door for tooling, automation and
   * inspectors (the u2 counterpart of the platform's `Widget.find`). Best-effort: the last
   * control constructed on (or assigned) a root wins it. */
  static forElement(el: Element | null): Control | undefined {
    for (let node = el; node; node = node.parentElement) {
      const control = Control._byRoot.get(node);
      if (control)
        return control;
    }
    return undefined;
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
