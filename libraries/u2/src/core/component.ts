import {Scope} from './scope.js';
import {Signal, isWritableSignal} from './signals.js';
import type {PropertyLike} from './property-like.js';
import type {BindPropLike, BindSourceLike, ComponentMetaLike, FuncLike, ObservableLike, WidgetLike,
  WidgetStatusLike} from './widget-like.js';

/** Base of every u2 component, visual or not: a name, an effect/cleanup scope, and the widget
 * introspection surface (DD7) generated from {@link meta}. A component constructed inside
 * another's scope is disposed with it. Platform-free; in Datagrok, `u2/dg`'s `host()` wires
 * {@link dispose} into the existing `DG.Widget` kill channel and delegates the introspection to a
 * real `DG.Widget`. Standalone hosts (gallery, tests) call {@link dispose} directly. */
export class Component implements WidgetLike, BindSourceLike {
  readonly scope: Scope;
  /** Unique within a spec: the anchor for selection, patches and automation. */
  name?: string;
  /** The registry metadata of the tag that built this component, stamped by the spec renderer —
   * the single source the introspection surface below is generated from. */
  meta?: ComponentMetaLike;
  /** The props the registry's `create` was handed, stamped alongside {@link meta}: the last resort
   * of the property read below, for a component that keeps a prop under another name or not at all
   * (a `create` wrapping a bare element). Construction-time values — a live one always wins. */
  specProps?: Record<string, unknown>;

  private _aiDescription: string | null = null;
  private _properties: PropertyLike[] | null = null;
  private _propertiesMeta?: ComponentMetaLike;
  private readonly _functions: FuncLike[] = [];
  private readonly _listeners: {id: string | null, next: (x: unknown) => void}[] = [];

  constructor(scope?: Scope) {
    this.scope = scope ?? new Scope();
    const owner = Scope.ambient;
    if (owner && owner !== this.scope)
      owner.own(() => this.dispose());
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
    return this._aiDescription ?? this.meta?.usage ?? this.meta?.description ?? null;
  }

  set aiDescription(x: string | null) {
    this._aiDescription = x;
  }

  /** The props {@link meta} declares, each closed over the component's own state: a same-named
   * signal reads and writes live, a same-named plain member, an `options` entry or — last —
   * a {@link specProps} entry reads read-only. {@link name} — the component's spec identity —
   * is never a property; an input's
   * declared `name` prop is its form key, a different thing that IS a property like any other.
   * A component never hand-writes its property list. The list is generated once per {@link meta}:
   * the accessors read live state, so only a re-stamp invalidates it. */
  getProperties(): PropertyLike[] {
    if (this._properties && this._propertiesMeta === this.meta)
      return this._properties;
    const self = this as unknown as Record<string, unknown>;
    this._propertiesMeta = this.meta;
    this._properties = (this.meta?.props ?? []).map((p) => {
      const prop: PropertyLike = {...p, get: () => Component._read(self, p.name)};
      if (self[p.name] instanceof Signal) {
        prop.set = (_source: unknown, value: unknown) => {
          (self[p.name] as Signal<unknown>).value = value;
        };
      }
      return prop;
    });
    return this._properties;
  }

  getFunctions(): FuncLike[] {
    return this._functions;
  }

  /** Events fired by this component, and — on a {@link Control} — the DOM events of that name
   * raised on its root. Without an id, every component event. The two streams are merged as they
   * are, so a subscriber sees either payload and must discriminate (see {@link ObservableLike}). */
  onEvent(eventId: string | null = null): ObservableLike<unknown> {
    return {
      subscribe: (next: (x: unknown) => void) => {
        const listener = {id: eventId, next};
        this._listeners.push(listener);
        const off = () => {
          const i = this._listeners.indexOf(listener);
          if (i >= 0)
            this._listeners.splice(i, 1);
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

  /** One binding step (DD5): a meta-declared prop's same-named Signal member — a named node is a
   * bind source, generated off {@link meta} the way {@link getProperties} is. '' answers the
   * default binding: the component's own `value` signal, where one exists. */
  bindStep(name: string): Signal<unknown> | BindSourceLike | null {
    const self = this as unknown as Record<string, unknown>;
    const member = name === '' ? self.value :
      this.meta?.props.some((p) => p.name === name) ? self[name] : undefined;
    return member instanceof Signal ? member as Signal<unknown> : null;
  }

  /** The signal-backed subset of the meta props — what a binding picker offers on a named node. */
  bindProps(): BindPropLike[] {
    const self = this as unknown as Record<string, unknown>;
    return (this.meta?.props ?? [])
      .filter((p) => self[p.name] instanceof Signal)
      .map((p) => ({...p, writable: isWritableSignal(self[p.name])}));
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
    this._functions.push(f);
  }

  protected fireEvent(eventId: string, args?: unknown): void {
    for (const listener of this._listeners.slice()) {
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
      root.append(child instanceof Control ? child.root : child);
    return new Control(root, scope);
  }
}
