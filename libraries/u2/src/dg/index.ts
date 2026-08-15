/* The ONLY layer allowed to import datagrok-api (enforced by eslint). Bridges u2's
   platform-free core onto the platform's existing lifecycle and event surfaces. */
export {asDartInput, DartInput} from './input-bridge.js';
export {columnInput, tableInput} from './pickers.js';
export {ObjectForm, propertyForm, objectForm} from './object-form.js';
export type {PropertyLike, PropertySource, ObjectFormOptions, FieldOverride} from './object-form.js';
export {fromDartInput, PlatformInput} from './from-dart-input.js';
export type {DartInputLike} from './from-dart-input.js';
export {userInput} from './user-input.js';
export {dapiSource, dapiPager, sanitizeFilterValue} from './dapi-source.js';
export type {DapiSourceLike, DapiSourceOptions, DapiPagerSourceLike, DapiPagerOptions} from './dapi-source.js';
export {handlerRenderer, HandlerRenderer, chip, EntityChip, entityCard, EntityCard, entityInput}
  from './entity.js';
export type {ChipOptions, EntityInputOptions} from './entity.js';
export {Editors} from './editors.js';
export type {EditorRule} from './editors.js';
export {moleculeInput, moleculeRenderer} from './molecule.js';
export type {MoleculeRendererOptions} from './molecule.js';

import * as DG from 'datagrok-api/dg';
import {Observable, Subscription} from 'rxjs';
import {Component} from '../core/component.js';
import {Scope} from '../core/scope.js';
import {signal, Signal, ReadonlySignal, rawEffect} from '../core/signals.js';

/** Hosts a u2 component as a `DG.Widget`: the component's disposal joins the widget's
 * `subs`, and `toDart()` registers it in the platform widget registry, so the Dart
 * kill-walk (view close) disposes all effects with zero new machinery.
 * Pass `closeIn` when docking the component ad-hoc: the pane ✕ fires only
 * `dockManager.onClosed` (never `Widget.kill` — platform gap, see PLAN.md P3), so the
 * auto-wire detaches on close. Components never actually mounted stay the caller's
 * responsibility to dispose. */
export function host(component: Component, closeIn?: DG.DockManager): DG.Widget {
  const w = DG.Widget.fromRoot(component.root);
  // core contract (2026-08): a docked pane's ✕ kills elements carrying this attribute;
  // inert on older cores, where the closeIn onClosed auto-wire below covers it
  component.root.setAttribute('data-kill-on-close', 'true');
  w.subs.push(new Subscription(() => component.dispose()));
  if (closeIn) {
    w.subs.push(closeIn.onClosed.subscribe((el: HTMLElement) => {
      if (el != null && (el === component.root || el.contains(component.root)))
        w.detach();
    }));
  }
  w.toDart();
  return w;
}

type ChromeItem = Component | HTMLElement;

export interface AppViewOptions {
  name: string;
  /** The view content — hosted as a `DG.Widget`, so view close disposes it (see {@link host}). */
  content: Component;
  /** Ribbon panel groups → `view.setRibbonPanels`. Main view controls (filter, refresh, mode
   * switches, menu bar) belong here, not inside the content area. */
  ribbon?: ChromeItem[][];
  /** → the shell's per-view status bar. A signal renders as a live text panel; elements and
   * components are placed as given. */
  status?: ReadonlySignal<string> | ChromeItem | ChromeItem[];
}

/** A platform view over a u2 component tree that rides the shell's chrome instead of
 * recreating it: ribbon for view controls, the per-view status bar for state. Chrome lives
 * outside the view root (the shell owns those containers), so the Dart kill-walk never
 * reaches it — components passed as chrome are disposed with the content instead. */
export function appView(options: AppViewOptions): DG.ViewBase {
  const w = host(options.content);
  const view = DG.View.fromRoot(options.content.root);
  view.name = options.name;

  const adopt = (item: ChromeItem): HTMLElement => {
    if (!(item instanceof Component))
      return item;
    w.subs.push(new Subscription(() => item.dispose()));
    return item.root;
  };
  const panel = (item: ChromeItem): HTMLDivElement => {
    const el = adopt(item);
    if (el instanceof HTMLDivElement)
      return el;
    const wrap = document.createElement('div');
    wrap.append(el);
    return wrap;
  };

  if (options.ribbon)
    view.setRibbonPanels(options.ribbon.map((group) => group.map(adopt)));
  const status = options.status;
  if (status instanceof Signal) {
    const el = document.createElement('div');
    w.subs.push(new Subscription(rawEffect(() => {
      el.textContent = status.value;
    })));
    view.statusBarPanels = [el];
  } else if (status != null) {
    const items = (Array.isArray(status) ? status : [status]) as ChromeItem[];
    view.statusBarPanels = items.map(panel);
  }
  return view;
}

/** Mirrors an rxjs observable into a signal owned by `scope`. */
export function toSignal<T>(observable: Observable<T>, initial: T, scope: Scope): ReadonlySignal<T> {
  const s = signal(initial);
  const sub = observable.subscribe((v) => s.value = v);
  scope.own(() => sub.unsubscribe());
  return s;
}

/** Exposes a signal as an rxjs observable (emits current value on subscribe).
 * The dual event surface: new code binds signals, old code subscribes — no migration. */
export function toObservable<T>(source: ReadonlySignal<T>): Observable<T> {
  return new Observable<T>((observer) => rawEffect(() => observer.next(source.value)));
}

/** Live u2 scopes vs platform-registered widgets — the leak detector's report.
 * After a view with u2 content closes, both numbers must return to their baseline. */
export function leakReport(): {liveScopes: number, liveWidgets: number} {
  return {liveScopes: Scope.liveCount, liveWidgets: DG.Widget.getAll().length};
}
