/** An input annotated `{editor: Package:func}` is edited with what that function returns, a `DG.InputBase`.
 *  The annotation is Flow's: the platform's run dialog keeps the plain input for the parameter's type.
 *  The editor function's inputs are matched by name to the node's inputs and get their values, tables
 *  included; those are the captured ones (a panel render never runs the flow), so until the flow has run
 *  up to the node the panel offers to run it. A list parameter stores the editor's value, a string list;
 *  any other type its `stringValue`. */

import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import type {CustomEditorContext, CustomInputEditor, CustomInputEditorFactory} from '../../utils/func-input-overrides';

/** The editor function an input names, or null. */
export function editorAnnotation(param: DG.Property): string | null {
  try {
    const ref = param.options?.['editor'];
    return typeof ref === 'string' && ref.trim() !== '' ? ref.trim() : null;
  } catch {
    return null; // Dart proxy reads can throw
  }
}

function isListParam(param: DG.Property): boolean {
  const t = String(param.propertyType);
  return t === 'list' || t === 'string_list';
}

export function annotatedEditor(ref: string): CustomInputEditorFactory {
  return (param: DG.Property, ctx: CustomEditorContext): CustomInputEditor => {
    const host = ui.div([], 'ff-annotated-editor');
    const list = isListParam(param);
    let current: unknown = null;
    let input: DG.InputBase | null = null;
    let subscription: {unsubscribe: () => void} | null = null;
    let detached = false;

    const valueOf = (i: DG.InputBase): unknown => list ? [...(i.value as unknown[] ?? [])] : i.stringValue;
    const show = (i: DG.InputBase, v: unknown): void => {
      if (v === null || v === undefined || v === '') return;
      if (list) i.value = Array.isArray(v) ? v : String(v).split(',').map((s) => s.trim()).filter(Boolean);
      else i.stringValue = String(v);
    };

    const ed: CustomInputEditor = {
      element: host,
      getValue: () => input ? valueOf(input) : current,
      setValue: (v) => {
        current = v;
        if (input) show(input, v);
      },
      isValid: () => input == null || input.validate(),
      detach: () => {
        detached = true;
        subscription?.unsubscribe();
      },
    };

    const notice = (text: string, action?: HTMLElement): void => {
      host.setAttribute('data-blocked', 'true');
      host.replaceChildren(ui.divV([
        ui.divH([ui.iconFA('table'), ui.divText(text, 'ff-annotated-editor-note')]),
        ...(action ? [action] : []),
      ], 'ff-annotated-editor-blocked'));
    };

    let func: DG.Func | null = null;
    try {
      func = DG.Func.byName(ref);
    } catch {/* not registered */}
    if (!func) {
      notice(`The editor ${ref} is not registered: is its package installed?`);
      return ed;
    }
    const editorFunc = func;
    const tables = editorFunc.inputs.filter((p) => String(p.propertyType) === 'dataframe').map((p) => p.name);
    const missing = (): string[] => tables.filter((t) => ctx.table?.(t) == null);

    const build = async (): Promise<void> => {
      const args: Record<string, unknown> = {};
      for (const p of editorFunc.inputs)
        args[p.name] = String(p.propertyType) === 'dataframe' ? ctx.table!(p.name) : ctx.inputValue(p.name);
      host.removeAttribute('data-blocked');
      host.replaceChildren(ui.loader());
      try {
        const made = await editorFunc.apply(args);
        if (detached) return;
        if (!(made instanceof DG.InputBase))
          throw new Error(`${ref} returned ${made == null ? 'nothing' : typeof made}, not an input`);
        input = made;
        show(input, current);
        subscription = input.onChanged.subscribe(() => ed.onChanged?.(valueOf(input!)));
        host.replaceChildren(input.root);
      } catch (e) {
        notice(`The editor failed: ${(e as Error)?.message ?? e}`);
      }
    };

    const unconnected = tables.filter((t) => ctx.isConnected && !ctx.isConnected(t));
    if (unconnected.length > 0)
      notice(`Connect ${unconnected.join(' and ')} to pick from ${unconnected.length === 1 ? 'it' : 'them'}.`);
    else if (missing().length === 0)
      void build();
    else {
      const run = ui.button('Run up to here', async () => {
        run.disabled = true;
        for (const t of missing()) await ctx.produceTable?.(t);
        if (detached) return;
        if (missing().length === 0) await build();
        else notice('The flow did not produce the tables this editor needs; run it and look at the errors.');
      });
      notice(`Picks from the ${tables.join(' and ')} this node receives, once the flow has run up to here.`, run);
    }
    return ed;
  };
}
