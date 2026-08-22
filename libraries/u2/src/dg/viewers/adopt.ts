/* The interim shape of "DG.Widget is a Control" (VIEWERS_AND_WIDGETS V2): the Control facets are
   mixed into the cached js-api wrapper's prototype chain, so the widget keeps its identity —
   `instanceof DG.Grid` holds, `DG.toJs(v.dart) === v`, events keep arriving on the same object —
   and every member the platform already has wins over the u2 one of that name. */
import * as DG from 'datagrok-api/dg';
import {Component, Control} from '../../core/component.js';
import type {ComponentState, PropertyChange} from '../../core/component.js';
import {signal} from '../../core/signals.js';
import type {BindPropLike, BindSourceLike, ObservableLike} from '../../core/widget-like.js';
import {dfBindings} from '../../sources/df-bindings.js';
import type {DataFrameLike} from '../../sources/df-bindings.js';
import {specWarn} from '../../spec/spec.js';

/** The kill-walk interop, reached through feature-detected globals read at call time: the running
 * js-api bundle predates the typed statics (`DG.Widget.registerCleanup`), and js-api edits are out
 * of scope. */
const api = globalThis as {
  grok_Widget_Kill?: (element: Element) => void;
  grok_Widget_RegisterCleanup?: (element: Element, cleanup: () => void) => void;
};

type Adopted = DG.Widget & Control & {_u2: ComponentState};
type Facets = Record<string, unknown> & ThisType<Adopted>;

const MIXED = new WeakMap<object, object>();
/** What a JS-owned widget's `onPropertyChanged` tells the property tier. */
const CHANGES = new WeakMap<object, Set<(change: PropertyChange) => void>>();

/** The walkable `table` step every adopted viewer answers first (VP-8): the DataFrame binding
 * surface over the viewer's frame, repointing with it. */
export const TABLE_STEP: BindPropLike = Object.freeze({name: 'table', type: 'dataframe', walkable: true});

/** Everything whose lifecycle is Dart's — every core viewer, every Dart widget — is released through
 * the kill-walk; a JS-implemented widget detaches both ways. */
export function dartOwned(widget: DG.Widget): boolean {
  return widget instanceof DG.DartWidget || (widget instanceof DG.Viewer && !(widget instanceof DG.JsViewer));
}

/** Makes `widget` a u2 {@link Control} in place and answers the same object: idempotent, one scope
 * per widget, adopted by the ambient scope like any component built inside one. A Dart-owned
 * widget is killed when its scope is disposed and its scope is disposed when the platform kills it;
 * a JS-owned one is detached when its scope is disposed and vice versa. A `DartWidget` announces no
 * property changes — its properties bind one-way until the platform offers a notification. */
export function adopt<T extends DG.Widget>(widget: T): T & Control {
  if (Control.is(widget))
    return widget;
  const proto = Object.getPrototypeOf(widget) as object;
  let mixed = MIXED.get(proto);
  if (mixed === undefined) {
    mixed = Object.create(proto, facets(proto, widget)) as object;
    MIXED.set(proto, mixed);
  }
  Object.setPrototypeOf(widget, mixed);
  const adopted = widget as unknown as Adopted;
  Component.init(adopted);
  adopted.root.classList.add('u2-component');
  if (widget instanceof DG.Viewer)
    adopted.root.classList.add('u2-viewer');
  if (dartOwned(widget))
    dartLifecycle(adopted);
  else
    jsLifecycle(adopted);
  return adopted as unknown as T & Control;
}

function descriptorOn(proto: object | null, name: string): PropertyDescriptor | undefined {
  for (let at = proto; at !== null; at = Object.getPrototypeOf(at)) {
    const found = Object.getOwnPropertyDescriptor(at, name);
    if (found)
      return found;
  }
  return undefined;
}

/** The u2 members the platform class lacks, plus the forced overrides of VP-6/VP-7 for its kind. */
function facets(proto: object, widget: DG.Widget): PropertyDescriptorMap {
  const out: PropertyDescriptorMap = {};
  for (const source of [Component.prototype, Control.prototype]) {
    for (const name of Object.getOwnPropertyNames(source)) {
      if (name !== 'constructor' && !(name in proto))
        out[name] = Object.getOwnPropertyDescriptor(source, name)!;
    }
  }
  const viewer = widget instanceof DG.Viewer;
  const dart = dartOwned(widget);
  const forced: Facets[] = [common(proto, viewer)];
  if (dart)
    forced.push(DART);
  if (viewer)
    forced.push(VIEWER);
  if (viewer && dart)
    forced.push(VIEWER_CHANGES);
  if (!dart)
    forced.push(jsOwned(proto as DG.Widget));
  for (const group of forced) {
    for (const name of Object.getOwnPropertyNames(group))
      out[name] = {...Object.getOwnPropertyDescriptor(group, name)!, enumerable: false};
  }
  return out;
}

function common(proto: object, viewer: boolean): Facets {
  const platformMeta = descriptorOn(proto, 'meta');
  return {
    get propertyTier() {
      return true;
    },
    // the registry stamps `component.meta = …`; a viewer's own `meta` is a getter-only helper (P8)
    get meta() {
      return this._u2.meta ?? platformMeta?.get?.call(this);
    },
    set meta(x) {
      this._u2.meta = x as ComponentState['meta'];
    },
    dispose() {
      this.scope.dispose();
    },
    bindStep(name: string) {
      if (viewer && name === 'table')
        return tableOf(this as Adopted & DG.Viewer);
      return Component.prototype.bindStep.call(this, name);
    },
  };
}

/** A viewer's property lives on its look, which the `props` bag is over (P6); a Dart widget's
 * takes the dart handle. */
const DART: Facets = {
  readProperty(name: string) {
    return this instanceof DG.Viewer ? (this.props as Record<string, unknown>)[name] :
      this.getProperties().find((p) => p.name === name)?.get?.(this.dart);
  },
  writeProperty(name: string, value: unknown) {
    if (this instanceof DG.FilterGroup && name === 'filters') {
      // the only `filters` write that rebuilds the panel (P14)
      for (const state of value as DG.FilterState[])
        this.updateOrAdd(state);
      return;
    }
    if (this instanceof DG.Viewer)
      (this.props as Record<string, unknown>)[name] = value;
    else
      this.getProperties().find((p) => p.name === name)?.set?.(this.dart, value);
  },
};

const VIEWER: Facets = {
  bindProps() {
    return [TABLE_STEP, ...Component.prototype.bindProps.call(this)];
  },
};

const VIEWER_CHANGES: Facets = {
  propertyChanges(): ObservableLike<PropertyChange> {
    const viewer = this as unknown as DG.Viewer;
    return {
      subscribe: (next) => {
        const subscription = viewer.onPropertyValueChanged.subscribe((e: unknown) =>
          next((e as {args?: {property?: {name?: string}}})?.args?.property?.name ?? null));
        return {unsubscribe: () => subscription.unsubscribe()};
      },
    };
  },
};

function jsOwned(proto: DG.Widget): Facets {
  return {
    propertyChanges(): ObservableLike<PropertyChange> {
      return {
        subscribe: (next) => {
          let listeners = CHANGES.get(this);
          if (listeners === undefined) {
            listeners = new Set();
            CHANGES.set(this, listeners);
          }
          listeners.add(next);
          return {unsubscribe: () => listeners!.delete(next)};
        },
      };
    },
    onPropertyChanged(property: DG.Property | null) {
      for (const next of CHANGES.get(this) ?? [])
        next(property?.name ?? null);
      proto.onPropertyChanged.call(this, property);
    },
    detach() {
      this._u2.detaching = true;
      this.scope.dispose();
      proto.detach.call(this);
    },
  };
}

function frameOf(viewer: DG.Viewer): DataFrameLike | undefined {
  return (viewer.dataFrame ?? undefined) as unknown as DataFrameLike | undefined;
}

function tableOf(viewer: Adopted & DG.Viewer): BindSourceLike {
  const state = viewer._u2;
  if (state.table === undefined) {
    const frame = signal(frameOf(viewer));
    const subscription = viewer.onDataFrameChanged.subscribe(() => frame.value = frameOf(viewer));
    viewer.scope.own(() => subscription.unsubscribe());
    state.table = dfBindings(frame, viewer.scope);
  }
  return state.table as BindSourceLike;
}

function dartLifecycle(widget: Adopted): void {
  if (api.grok_Widget_Kill === undefined || api.grok_Widget_RegisterCleanup === undefined) {
    specWarn('viewers: grok_Widget_Kill / grok_Widget_RegisterCleanup are not available on this ' +
      'client — an adopted Dart widget is not released with its scope');
    return;
  }
  const state = widget._u2;
  const root = widget.root;
  api.grok_Widget_RegisterCleanup(root, () => {
    state.platformKilled = true;
    widget.scope.dispose();
  });
  widget.scope.own(() => {
    if (!state.platformKilled)
      api.grok_Widget_Kill?.(root);
  });
}

function jsLifecycle(widget: Adopted): void {
  widget.scope.own(() => {
    if (!widget.isDetached && !widget._u2.detaching)
      widget.detach();
  });
}
