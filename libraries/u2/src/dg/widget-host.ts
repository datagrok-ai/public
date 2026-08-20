/* The platform boundary of the widget convergence (DD7): a hosted u2 component IS a `DG.Widget`,
   whose whole introspection surface delegates to it. The assertion below keeps the core's
   structural `WidgetLike` from drifting off the real interface it copies. */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {Observable, Subscription} from 'rxjs';
import {Control} from '../core/component.js';
import type {PropertyLike} from '../core/property-like.js';
import type {FuncLike, WidgetLike} from '../core/widget-like.js';

const _conformance: WidgetLike = null! as unknown as DG.Widget;

/** Namespace the ad-hoc component functions are minted into, so they never land in `core` and
 * override a platform function of the same name. */
const FUNC_NAMESPACE = 'U2';

/** Registry-only prop flags (`SpecPropMeta`): a platform property has nowhere to put them. */
const U2_ONLY_FLAGS = new Set(['bindable', 'twoWay']);

/** A u2 component as the platform sees it: `getProperties`, `getFunctions`, `onEvent`,
 * `aiDescription` and `getWidgetStatus` are the component's own, and the component is disposed
 * through the widget's kill channel. */
export class U2Widget extends DG.Widget {
  readonly component: Control;

  constructor(component: Control) {
    super(component.root);
    this.component = component;
  }

  getProperties(): DG.Property[] {
    if (this._properties.length === 0)
      this._properties = this.component.getProperties().map((p) => U2Widget._property(p));
    return this._properties;
  }

  /** The component's callables as platform functions, minted on first ask and cached. The mint is
   * a global registration, so a component that declares none registers nothing. */
  getFunctions(): DG.Func[] {
    const declared = this.component.getFunctions();
    if (this._functions.length === 0 && declared.length > 0)
      this._functions = declared.map((f) => U2Widget._func(this.component, f));
    return this._functions;
  }

  get aiDescription(): string | null {
    return this.component.aiDescription;
  }

  set aiDescription(x: string | null) {
    this.component.aiDescription = x;
  }

  onEvent(eventId: string | null = null): Observable<any> {
    return new Observable<any>((observer) => {
      const subscription = this.component.onEvent(eventId).subscribe((x) => observer.next(x));
      return () => subscription.unsubscribe();
    });
  }

  getWidgetStatus(): DG.IWidgetStatus {
    return this.component.getWidgetStatus() as DG.IWidgetStatus;
  }

  private static _property(p: PropertyLike): DG.Property {
    const {name, type, propertyType, get, set, ...rest} = p;
    const options: Record<string, unknown> = {};
    for (const [key, value] of Object.entries(rest)) {
      if (value != null && !U2_ONLY_FLAGS.has(key))
        options[key] = value;
    }
    const property = DG.Property.create(name, (propertyType ?? type ?? DG.TYPE.OBJECT) as DG.Type,
      () => get?.(null), (set ?? null) as unknown as DG.PropertySetter, null);
    return property.fromOptions(options);
  }

  /** `grok.functions.register` is the only client-side door to a real `Func` (js-api
   * `functions.ts:137`, the `DomainGrid.widgetActionFunc` precedent) — it takes a signature and a
   * closure, and persists nothing. Params arrive positionally, in the order declared.
   * Known retention: there is no unregister, so a minted Func pins its component for the session,
   * and a second registration of the same name rebinds name lookups to the newer component while
   * each widget keeps its own handle. The component's `name` goes into the minted name to keep
   * distinct instances apart. */
  private static _func(component: Control, f: FuncLike): DG.Func {
    const owner = component.meta?.tag ?? component.constructor.name;
    const scope = component.name ? `${owner}_${component.name}` : owner;
    const params = f.inputs.map((p) => `${p.propertyType ?? p.type ?? DG.TYPE.OBJECT} ${p.name}`);
    return grok.functions.register({
      signature: `${DG.TYPE.OBJECT} ${U2Widget._identifier(scope)}_${U2Widget._identifier(f.name)}` +
        `(${params.join(', ')})`,
      run: async (...args: unknown[]) => {
        const values: Record<string, unknown> = {};
        for (let i = 0; i < f.inputs.length; i++)
          values[f.inputs[i].name] = args[i];
        return await f.apply(values);
      },
      isAsync: true,
      namespace: FUNC_NAMESPACE,
      options: f.description ? {description: f.description} : {},
    });
  }

  private static _identifier(name: string): string {
    return name.replace(/[^A-Za-z0-9]/g, '_');
  }
}

/** Hosts a u2 component as a `DG.Widget`: the component's disposal joins the widget's
 * `subs`, and `toDart()` registers it in the platform widget registry, so the Dart
 * kill-walk (view close) disposes all effects with zero new machinery.
 * Pass `closeIn` when docking the component ad-hoc: the pane ✕ fires only
 * `dockManager.onClosed` (never `Widget.kill` — platform gap, see PLAN.md P3), so the
 * auto-wire detaches on close. Components never actually mounted stay the caller's
 * responsibility to dispose. */
export function host(component: Control, closeIn?: DG.DockManager): DG.Widget {
  const w = new U2Widget(component);
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
