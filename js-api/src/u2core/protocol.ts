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

/** JSON-safe projection of a value — what crosses `page.evaluate`. Dates become ISO strings,
 * entities/objects with a `name` become that name, arrays map element-wise, plain objects key-wise,
 * functions/elements become undefined, everything else `String(x)`. Produced by `Control.jsonValue`. */
export type JsonValue = null | boolean | number | string | JsonValue[] | {[key: string]: JsonValue};

/** What the user does to a control or an item — the bdd library's gesture verbs. */
export type U2Gesture =
  | 'type' | 'clear' | 'pick' | 'toggle' | 'press' | 'expand' | 'collapse' | 'select'
  | 'open' | 'close' | 'drag' | 'resize' | 'hover' | 'scroll' | 'add' | 'remove' | 'rename';

/** One addressable item of a collection-like control: a list/tree/table row, a chip, an option,
 * a radio/checkbox item, a tab, a wizard step, a menu item, a breadcrumb. */
export interface IItemState {
  /** Stable identity, equal to the item element's `data-u2-item`. */
  key: string;
  /** The text the user reads on the item (trimmed, single-spaced). */
  label: string;
  /** Position among ALL items (not the rendered window), 0-based. */
  index: number;
  /** Selectable items: in the selection (`aria-selected`). */
  selected?: boolean;
  /** Radio / checkbox / checkable tree items (`aria-checked`). */
  checked?: boolean | 'mixed';
  /** Tree nodes (`aria-expanded`). */
  expanded?: boolean;
  /** Tree depth, 1-based (aria-level). */
  level?: number;
  /** The item cannot be picked (`aria-disabled`). */
  disabled?: boolean;
  /** The current tab / the current wizard step. */
  current?: boolean;
  /** Wizard step already passed. */
  completed?: boolean;
  /** Wizard step / row reporting an error. */
  invalid?: boolean;
  /** Whether a DOM element exists for it right now (false: off-screen in a virtual list). */
  rendered: boolean;
}

/** The items of a collection-like control. */
export interface IItemsState {
  /** All items, not the rendered window. */
  count: number;
  /** The first {@link ITEMS_CAP} items; `truncated` says whether there are more. */
  list: IItemState[];
  /** More than {@link ITEMS_CAP} items, and `list` holds the first ones only. */
  truncated: boolean;
}

/** The typed, JSON-safe state block every `Control` fills — what generic bdd steps read.
 * Keys a kind does not have are absent (never null-filled). Read on demand, never cached. */
export interface IControlState {
  /** The `data-u2` value; for a platform `DG.Widget` without one, `widget.type` ('Scatter plot'). */
  kind: string;
  /** The registry tag when the control was built from a spec node or its kind maps to one. */
  tag?: string;
  /** Explicit identity — equal to `data-u2-name`. Never derived from a caption. */
  name?: string;
  /** The accessible name the user sees: label part, button text, dialog/section title. */
  label?: string;
  /** No `aria-disabled="true"` / `disabled` on the root. */
  enabled: boolean;
  /** Why it is disabled — the root's `aria-description` while it is; absent when enabled or disabled without a
   * reason. */
  disabledReason?: string;
  /** Connected, not under a `hidden` element, and laid out. */
  visible: boolean;
  /** `document.activeElement` is the root or inside it (or inside a popup this control owns). */
  focused: boolean;
  /** Shown but not editable (`aria-readonly`). */
  readOnly?: boolean;
  /** Must have a value (`aria-required`). */
  required?: boolean;
  /** Inputs, forms, wizard steps: no validation message. */
  valid?: boolean;
  /** The validation / error-row message the user sees; null when none. */
  error?: string | null;
  /** Async work in flight (dialog OK, async view, source load, wizard step). */
  busy?: boolean;
  /** Raw value, JSON-projected (inputs; a range slider: `[lo, hi]`). */
  value?: JsonValue;
  /** The value as displayed (the editor's text, the chosen option's caption, a checkbox 'true'/'false'). */
  text?: string;
  /** The hint the empty editor shows. */
  placeholder?: string;
  /** Popup / dropdown / dialog open. */
  open?: boolean;
  /** Section / accordion pane / tree node / disclosure expanded. */
  expanded?: boolean;
  /** Checkbox / switch state. */
  checked?: boolean | 'mixed';
  /** Toggle button state (`aria-pressed`). */
  pressed?: boolean;
  /** For a selectable card / row-like control. */
  selected?: boolean;
  /** Collection-like controls: every item, rendered or not. */
  items?: IItemsState;
  /** Keys of the selected items (list, tree, table, multi-choice, chips, multi-select). */
  selection?: string[];
  /** Key of the current tab / wizard step / active breadcrumb. */
  current?: string;
  /** Free per-kind readings that do not deserve a typed key (a stat card's delta). Same contract as
   * {@link IWidgetStatus.values}: number | string | boolean only. */
  extra?: {[key: string]: number | string | boolean};
}

/** Runtime snapshot of a widget's structure, read by the automated testing system. */
export interface IWidgetStatus {
  parts: {[name: string]: Element};
  hitAreas: {[name: string]: IRectBounds};
  /** Named readings of the widget's state a test compares before and after a change (a color
   * scale's shown range). */
  values?: {[name: string]: number | string | boolean};
  shortcuts: {[key: string]: string};
  events: IEventType[];
  description: string | null;
  error: string | null;
  inputs?: IInputStatus[];
  /** Filled by every u2 `Control` and by `DG.Widget` (DOM-derived). Optional only because the Dart
   * status object does not carry it; every JS read path fills it. */
  state?: IControlState;
}

/** Item cap of {@link IItemsState.list}. */
export const ITEMS_CAP = 500;

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
