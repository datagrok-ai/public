/* The introspection and binding vocabulary of the u2 core — event streams, widget status,
   functions, the DD5 bind protocol and registry metadata. Structural on purpose: the core stays
   platform-free, and a platform `DG.Widget` answers all of it by inheritance. The surviving
   structural shapes each admit a wider population than any platform class — that is why they
   are not the class. */
import type {IProperty} from './property-fields.js';
import type {Signal} from './signals.js';

/** The subscription surface of an event stream; an rxjs `Observable` satisfies it structurally.
 * Deliberately NOT rxjs `Observable`: the notification seams return hand-rolled literals, and
 * the dependency-free core never types against rxjs (a webpack external).
 * What {@link Component.onEvent} emits is a union: a component event's own payload, or a DOM
 * `Event` where the control's root raises one of the same name. Nothing wraps or tags them —
 * a subscriber that cares tests for `instanceof Event`. */
export interface ObservableLike<T> {
  subscribe(next: (x: T) => void): {unsubscribe(): void};
}

/** Event type descriptor a widget reports. */
export interface IEventType {
  name: string;
  eventName: string;
  description: string;
}

/** Bounding rectangle of a named hit area. */
export interface IRectBounds {
  x: number;
  y: number;
  width: number;
  height: number;
}

/** Live state of ONE named input of a widget — the machine-readable counterpart of what the
 * user sees in a form. Values and validation are read on demand, never cached; `name` is the
 * name the widget's properties address the input by, so `props[name] = value` writes exactly
 * what typing into it would. */
export interface IInputStatus {
  name: string;
  caption?: string;
  type: string;
  semType?: string;
  value: any;
  choices?: any[];
  required: boolean;
  valid: boolean;
  error?: string;
  description?: string;
  ref?: string;
}

/** Runtime snapshot of a widget's structure, read by the automated testing system. */
export interface IWidgetStatus {
  parts: {[name: string]: Element};
  hitAreas: {[name: string]: IRectBounds};
  shortcuts: {[key: string]: string};
  events: IEventType[];
  description: string | null;
  error: string | null;
  inputs?: IInputStatus[];
}

/** A property whose name is guaranteed — what descriptors and bind steps carry. */
export interface NamedProperty extends IProperty {
  name: string;
}

/** A callable a component exposes, params as properties. The invocation is named after
 * `DG.Func.apply` so that a real platform function IS a `FuncLike` — the host hands out either.
 * Deliberately NOT `DG.Func`: components mint plain callables no registry knows. */
export interface FuncLike {
  name: string;
  description?: string;
  inputs: NamedProperty[];
  /** May return a promise. */
  apply(params?: Record<string, unknown>): unknown;
}

/** One enumerated binding step: a property (typed, semType included) plus what the picker
 * needs — no probing, no signal creation during enumeration. */
export interface BindProp extends NamedProperty {
  /** `bindStep(name)` answers another source — the picker renders an expandable node. */
  walkable?: boolean;
  /** The leaf is a writable signal — two-way capable. */
  writable?: boolean;
  /** What `bindStep('')` answers — a path that stops at the source binds here. */
  default?: boolean;
}

/** The DD5 binding protocol: a source resolves one path step at a time and enumerates what a
 * picker may offer without allocating a single signal — laziness is the contract, not an
 * optimization. Every {@link Component} answers it off its registry meta. */
export interface BindSource {
  /** Resolve one step; a Signal ends the walk, a source continues it, null is unresolvable.
   * '' asks for the default binding. */
  bindStep(name: string): Signal<unknown> | BindSource | null;
  /** What a binding picker enumerates. Must not allocate signals. */
  bindProps(): BindProp[];
}

/** The registry metadata a component introspects itself through; the spec layer's
 * `ComponentMeta` extends it (structural, so the core never imports the spec layer). */
export interface ComponentMetaBase {
  tag?: string;
  props: NamedProperty[];
  description?: string;
  usage?: string;
  events?: string[];
  /** The tag-level design-time sample a source falls back to (DD9). */
  designPreview?: object;
}
