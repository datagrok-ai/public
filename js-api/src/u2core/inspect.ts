/* The in-page accessors automation reads u2 controls through (`DG.U2.state(el)` …): resolve the
   control owning an element, answer its status, state block, items and declared properties. They
   resolve u2 controls only — a platform `DG.Widget` answers `DG.Widget.find(el)?.getWidgetStatus()`. */
import {Control} from './component.js';
import type {IControlState, IItemsState, IWidgetStatus, JsonValue} from './protocol.js';

/** One declared property of a control, as {@link props} reports it. */
export interface PropState {
  /** The property name `setProp` takes. */
  name: string;
  /** The property type (`string`, `bool`, `int` …); `object` when undeclared. */
  type: string;
  /** The current value, JSON-projected. */
  value: JsonValue | undefined;
  /** Whether `setProp` can write it. */
  writable: boolean;
}

/** The u2 control owning `el`: the nearest control root at or above it, or the owner of the popup it
 * sits in. Platform widget wrappers are passed over — whether one exists depends on who asked the
 * platform for it. Undefined for an element outside any u2 control.
 * @example
 * const input = DG.U2.control(document.querySelector('[data-u2="text-input"]'));
 * @see {@link https://public.datagrok.ai/js/samples/ui/widgets/u2-status} */
export function control(el: Element | null): Control | undefined {
  for (let found = Control.forElement(el); found !== undefined; found = Control.forElement(found.root.parentElement)) {
    if (typeof (found as {toDart?: unknown}).toDart !== 'function')
      return found;
  }
  return undefined;
}

/** Full status of the control owning `el`; its parts are elements, so use it inside the page only.
 * @example
 * const parts = DG.U2.status(el)?.parts;
 * @see {@link https://public.datagrok.ai/js/samples/ui/widgets/u2-status} */
export function status(el: Element | null): IWidgetStatus | null {
  return control(el)?.getWidgetStatus() ?? null;
}

/** JSON-safe state of the control owning `el` — what `page.evaluate` hands to a test runner; null
 * outside a u2 control. Synchronous and side-effect free.
 * @example
 * const {label, value, valid, enabled} = DG.U2.state(el)!;
 * @see {@link https://public.datagrok.ai/js/samples/ui/widgets/u2-status} */
export function state(el: Element | null): IControlState | null {
  return control(el)?.getState() ?? null;
}

/** Items of the control owning `el` (the state's `items`), uncapped with `all: true`; null when it
 * lists none.
 * @example
 * const keys = DG.U2.items(el, {all: true})?.list.map((item) => item.key);
 * @see {@link https://public.datagrok.ai/js/samples/ui/widgets/u2-status} */
export function items(el: Element | null, opts?: {all?: boolean}): IItemsState | null {
  return control(el)?.getItems(opts) ?? null;
}

/** Brings an item of the control owning `el` into view by key or label (see `Control.revealItem`);
 * rejects naming the candidates when a label matches several.
 * @example
 * const row = await DG.U2.revealItem(el, 'Discovery');
 * @see {@link https://public.datagrok.ai/js/samples/ui/widgets/u2-status} */
export async function revealItem(el: Element | null, keyOrLabel: string): Promise<HTMLElement | null> {
  return await control(el)?.revealItem(keyOrLabel) ?? null;
}

/** The declared properties of the control owning `el`, values JSON-projected.
 * @example
 * const writable = DG.U2.props(el).filter((p) => p.writable).map((p) => p.name);
 * @see {@link https://public.datagrok.ai/js/samples/ui/widgets/u2-status} */
export function props(el: Element | null): PropState[] {
  const c = control(el);
  if (c === undefined)
    return [];
  return c.getProperties().map((p) => ({name: p.name, type: p.type ?? p.propertyType ?? 'object',
    value: Control.jsonValue(p.get?.(c.propertyTarget)), writable: p.set != null}));
}

/** Writes a declared property of the control owning `el` through its setter — for putting a control in
 * a state from a test. Throws when the control has no such property or it is read-only.
 * @example
 * DG.U2.setProp(el, 'enabled', false);
 * @see {@link https://public.datagrok.ai/js/samples/ui/widgets/u2-status} */
export function setProp(el: Element | null, name: string, value: unknown): void {
  const c = control(el);
  if (c === undefined)
    throw new Error('u2: not inside a u2 control');
  const kind = c.root.dataset.u2 ?? c.componentMeta?.tag ?? 'control';
  const prop = c.getProperties().find((p) => p.name === name);
  if (prop === undefined)
    throw new Error(`u2: ${kind} has no property ${name}`);
  if (prop.set == null)
    throw new Error(`u2: ${kind}.${name} is read-only`);
  prop.set(c.propertyTarget, value);
}
