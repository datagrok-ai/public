/* The introspection surface every platform widget answers — properties, functions, events, an AI
   briefing and a status snapshot (DD7). Mirrored structurally, as `PropertyLike` and
   `DartInputLike` are, so the core stays platform-free while a hosted component IS a `DG.Widget`;
   src/dg/shell/widget-host.ts asserts at compile time that the real interface satisfies this one. */
import type {PropertyLike} from './property-like.js';
import type {Signal} from './signals.js';

/** The subscription surface of an event stream; an rxjs `Observable` satisfies it structurally.
 * What {@link WidgetLike.onEvent} emits is a union: a component event's own payload, or a DOM
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

/** A callable a component exposes, params as properties. The invocation is named after
 * `DG.Func.apply` so that a real platform function IS a `FuncLike` — the host hands out either. */
export interface FuncLike {
  name: string;
  description?: string;
  inputs: PropertyLike[];
  /** May return a promise. */
  apply(params?: Record<string, unknown>): unknown;
}

/** What the shell, the context panel, automation and the copilot ask of a widget — and, through
 * this interface, of any u2 component, visual or not. `DG.Widget` implements it as declared. */
export interface WidgetLike {
  getProperties(): PropertyLike[];
  getFunctions(): FuncLike[];
  onEvent(eventId?: string | null): ObservableLike<unknown>;
  aiDescription: string | null;
  getWidgetStatus(): IWidgetStatus;
}

/** Mirrors the spec layer's `BindProp` (structural, as {@link ComponentMetaLike} is): one
 * enumerated binding step, typed like a property plus what the picker needs. */
export interface BindPropLike extends PropertyLike {
  /** `bindStep(name)` answers another source — the picker renders an expandable node. */
  walkable?: boolean;
  /** The leaf is a writable signal — two-way capable. */
  writable?: boolean;
  /** What `bindStep('')` answers — a path that stops at the source binds here. */
  default?: boolean;
}

/** Mirrors the spec layer's `BindSource` — the DD5 binding protocol — structurally, so the core
 * never imports the spec layer. Every {@link Component} answers it off its registry meta. */
export interface BindSourceLike {
  /** Resolve one step; a Signal ends the walk, a source continues it, null is unresolvable.
   * '' asks for the default binding. */
  bindStep(name: string): Signal<unknown> | BindSourceLike | null;
  /** What a binding picker enumerates. Must not allocate signals. */
  bindProps(): BindPropLike[];
}

/** The registry metadata a component introspects itself through; `ComponentMeta` satisfies it
 * (structural, so the core never imports the spec layer). */
export interface ComponentMetaLike {
  tag?: string;
  props: PropertyLike[];
  description?: string;
  usage?: string;
  events?: string[];
  /** The tag-level design-time sample a source falls back to (DD9). */
  designPreview?: object;
}
