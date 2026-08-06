/**
 * The widget-action vocabulary: the domain widgets' actions as REAL platform
 * {@link DG.Func}s, so `widget.getFunctions()` returns something typed, annotated,
 * discoverable and invocable through the machinery the platform already has (and
 * already exposes to the AI integrations, `js-api/src/widgets/base.ts`).
 *
 * ONE shared vocabulary, not one function per widget instance: `Save`, `Discard`
 * and `Reset` each take the widget they act on (`//input: widget widget`) and call
 * the member of the same name on it, so a hundred open forms register nothing and
 * every one of them reports the same three names. Registration happens on first
 * use and is idempotent — the Func is cached per qualified name.
 *
 * ```ts
 * const save = form.getFunctions().find((f) => f.name === 'Save')!;
 * await save.apply({widget: form});
 * ```
 *
 * @module actions
 */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

/** Namespace every widget-action Func is registered under — what keeps `Save`
 * from colliding with another package's `Save`. */
export const WIDGET_ACTION_NAMESPACE = 'DomainUi';

/** Tags every widget-action Func carries: `domain-ui` marks the origin,
 * `widget-action` the calling convention (one `widget` input). */
export const WIDGET_ACTION_TAGS = 'domain-ui,widget-action';

/** What a widget must implement to be acted on by the vocabulary. Every member is
 * optional: a widget offers the actions it can perform (see
 * {@link DomainForm.getFunctions}), and one that cannot is simply not offered. */
export interface IWidgetActions {
  save?(): Promise<boolean>;
  discard?(): void | Promise<void>;
  reset?(): void | Promise<void>;
  addRow?(): number | Promise<number>;
  deleteRow?(): boolean | Promise<boolean>;
  refresh?(): void | Promise<any>;
}

const _funcs = new Map<string, DG.Func>();

/**
 * The widget-action Func called [name], registering it on first use.
 *
 * Idempotent by qualified name (`DomainUi:<name>`): the first registration wins and
 * every later call returns the same Func, so a library bundled into several
 * packages registers one vocabulary, not one per bundle.
 */
export function widgetActionFunc(name: string, run: (widget: any) => Promise<boolean>,
  description: string): DG.Func {
  const existing = _funcs.get(name);
  if (existing != null)
    return existing;
  const func = grok.functions.register({
    signature: `bool ${name}(widget widget)`,
    run: run,
    tags: WIDGET_ACTION_TAGS,
    isAsync: true,
    namespace: WIDGET_ACTION_NAMESPACE,
    options: {description: description},
  });
  _funcs.set(name, func);
  return func;
}

/** Writes the widget's pending changes (one transaction per editor); resolves to
 * whether everything landed. */
export function saveFunc(): DG.Func {
  return widgetActionFunc('Save', async (widget: any) =>
    typeof widget?.save === 'function' ? await widget.save() === true : false,
  'Saves the pending changes of the widget');
}

/** Drops the widget's pending changes. */
export function discardFunc(): DG.Func {
  return widgetActionFunc('Discard', async (widget: any) => {
    if (typeof widget?.discard !== 'function')
      return false;
    await widget.discard();
    return true;
  }, 'Discards the pending changes of the widget');
}

/** Puts the widget back to the values it was opened with. */
export function resetFunc(): DG.Func {
  return widgetActionFunc('Reset', async (widget: any) => {
    if (typeof widget?.reset !== 'function')
      return false;
    await widget.reset();
    return true;
  }, 'Restores the values the widget was opened with');
}

/** Appends an empty row to a grid's pending batch (prefilled with the grid's
 * `defaults`); false when the widget cannot insert. */
export function addRowFunc(): DG.Func {
  return widgetActionFunc('AddRow', async (widget: any) => {
    if (typeof widget?.addRow !== 'function')
      return false;
    return await widget.addRow() >= 0;
  }, 'Appends a new row to the widget');
}

/** Marks the widget's current row for deletion (the save writes it). */
export function deleteRowFunc(): DG.Func {
  return widgetActionFunc('DeleteRow', async (widget: any) => {
    if (typeof widget?.deleteRow !== 'function')
      return false;
    return await widget.deleteRow() === true;
  }, 'Marks the current row of the widget for deletion');
}

/** Re-reads the widget's rows from the server (through the unsaved-changes gate
 * of whatever owns it). */
export function refreshFunc(): DG.Func {
  return widgetActionFunc('Refresh', async (widget: any) => {
    if (typeof widget?.refresh !== 'function')
      return false;
    await widget.refresh();
    return true;
  }, 'Reloads the widget from the server');
}

/** A function name for [caption]: identifier characters only, so an action called
 * 'Edit...' registers as `Edit`. */
export function actionFuncName(caption: string): string {
  const name = `${caption ?? ''}`.replace(/[^A-Za-z0-9_]/g, '');
  return name === '' ? 'Action' : name;
}

/** Wraps a custom {@link DG.DomainAction} as a widget-action Func, through the same
 * idempotent registration — a plugin's extra action is as discoverable as `Save`. */
export function domainActionFunc(action: DG.DomainAction): DG.Func {
  return widgetActionFunc(actionFuncName(action.name), async () => {
    const result = await action.run();
    return result !== false;
  }, `${action.name}`);
}
