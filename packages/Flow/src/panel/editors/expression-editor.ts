/** Formula editor for the string parameters that hold an expression — the row
 *  conditions of Filter / Delete / Extract / Select Rows, and Expression To
 *  Column's formula.
 *
 *  A condition typed into a bare text field is guesswork: no column names, no
 *  function list, no validation until the node runs. The platform already
 *  solves this — the **Add New Column** editor — and core already reuses it as
 *  a *property* editor: every Dart viewer's `filter` property opens it through
 *  `PropertyViewFormulaEditor` (xamgle/property_grid/editors), flagged
 *  `aux.filterFormulaEditor` so it drops the name/type inputs, validates as
 *  `bool`, and never appends a column. This is the same idea, mounted inline
 *  and reactive instead of behind an ellipsis button.
 *
 *  **The table gates it.** The editor is built against a real DataFrame — it
 *  autocompletes that table's columns. Until one is available the parameter is
 *  an ordinary string input, with the reason spelled out and (when the upstream
 *  is wired but not yet computed) a button to run the slice. That mirrors the
 *  column picker's ladder; the difference is that the picker runs on a click
 *  while this renders on every panel render, so the table is only ever *read*
 *  from what has already been captured — never computed as a side effect.
 *
 *  **Reactivity** is PowerPack's: the widget publishes the debounced text onto
 *  the `AddNewColumn` call it was handed, and we read it back through
 *  `FuncCallParam.onChanged`.
 *
 *  **Validity** is PowerPack's too. Deciding whether a formula is a true/false
 *  condition means EVALUATING it — the platform infers a formula's type from
 *  its result — so the check is asynchronous and belongs where the computation
 *  already happens. The widget publishes each verdict on its own root (the only
 *  channel that survives the package crossing: a widget arrives here as a fresh
 *  JS wrapper around the same Dart handle); this module caches it per node so
 *  {@link expressionRequirements} can answer the node's readiness check
 *  synchronously. */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import type {
  CustomEditorContext, CustomInputEditor, CustomInputEditorFactory, FuncNodeValidator,
} from '../../utils/func-input-overrides';
import {getParamDisplayName} from '../../utils/dart-proxy-utils';
import {propertyNameToFriendly} from '../../utils/naming';

/** The PowerPack entry point that returns the formula editor as a widget.
 *  Absent (older PowerPack, or none deployed) → the plain string input, which
 *  is exactly the no-table state, so nothing needs a second fallback path. */
const EDITOR_FUNC = {package: 'PowerPack', name: 'expressionEditorWidget'};

/** PowerPack's validity contract (`add-new-column.ts`), duplicated as literals
 *  because packages don't import each other's source. */
const VALIDATED_EVENT = 'expression-validated';

/** The last verdict the editor received, per node — runtime-only and weakly
 *  held, deliberately NOT a node property: serializing it would dirty a flow
 *  nobody edited the moment a check lands. Same treatment as the MPO profile
 *  cache, for the same reasons. */
const validated = new WeakMap<object, {expression: string; error: string}>();

/** Readiness check for a parameter edited with this editor: the node is not
 *  runnable while its formula is known to be wrong.
 *
 *  Fails **OPEN** whenever it cannot decide — nothing cached, or the text has
 *  moved on since the last verdict. A verdict needs the upstream table, which
 *  only an open panel has; a flow loaded from disk and never opened would
 *  otherwise be permanently unrunnable. The op itself refuses a bad condition
 *  at run time, so an unknown state costs a clear error, not a wrong result. */
export function expressionRequirements(paramName: string): FuncNodeValidator {
  const label = propertyNameToFriendly(paramName);
  return (node) => {
    const verdict = validated.get(node);
    if (!verdict || verdict.error === '') return [];
    if (verdict.expression !== String(node.inputValues[paramName] ?? '')) return [];
    // First sentence only — the full text sits next to the editor already, and
    // this one has to fit on a node's hint line.
    return [`${label} — ${verdict.error.split('. ')[0]}`];
  };
}

export interface ExpressionEditorOptions {
  /** The dataframe parameter whose table the formula is written against. */
  tableParam: string;
  /** Constrain to a true/false formula (`aux.filterFormulaEditor`) — what makes
   *  it a *condition* editor rather than a general one. */
  booleanOnly: boolean;
}

/** Builds the factory for one parameter. Curried because
 *  `CUSTOM_FUNC_INPUT_EDITORS` maps a parameter to a bare factory, and the two
 *  flavours differ only in these options. */
export function expressionEditor(options: ExpressionEditorOptions): CustomInputEditorFactory {
  return (param, ctx) => buildExpressionEditor(param, ctx, options);
}

/** Condition editor: a boolean formula over the `table` input. */
export const rowConditionEditor = expressionEditor({tableParam: 'table', booleanOnly: true});

/** General formula editor: any expression over the `table` input. */
export const columnFormulaEditor = expressionEditor({tableParam: 'table', booleanOnly: false});

function buildExpressionEditor(
  param: DG.Property, ctx: CustomEditorContext, options: ExpressionEditorOptions,
): CustomInputEditor {
  const ed: CustomInputEditor = {} as CustomInputEditor;
  const host = ui.div([], 'ff-expression-editor');
  const label = getParamDisplayName(param);
  let value = '';
  /** The table the mounted formula editor was built against — rebuilding it on
   *  every render would throw away the user's cursor and undo history. */
  let builtFor: DG.DataFrame | null = null;
  let disposeWidget: (() => void) | null = null;

  const emit = (next: string): void => {
    if (next === value) return;
    value = next;
    ed.onChanged?.(next);
  };

  /** Store a verdict and make the rest of the UI act on it: the host element
   *  marks itself (styling and tests), and the node recomputes its hint and its
   *  run gate — nothing else re-renders a node because a check came back. */
  const recordVerdict = (detail: {expression: string; error: string}): void => {
    const error = String(detail?.error ?? '');
    host.toggleAttribute('data-invalid', error !== '');
    const node = ctx.node;
    if (!node) return;
    // Compare against "no verdict yet" as CLEAN, not as different-from-clean:
    // every mount validates, so treating the first clean answer as a change
    // reported a parameter edit on merely opening the panel — invalidating
    // downstream results, and (worse) rebuilding the panel out from under the
    // dialog this editor had just opened, leaving it writing to a dead call.
    const previousError = validated.get(node)?.error ?? '';
    validated.set(node, {expression: String(detail?.expression ?? ''), error});
    if (previousError !== error) node.editorBridge?.notifyParamsChanged(node.id);
  };

  const reset = (): void => {
    disposeWidget?.();
    disposeWidget = null;
    builtFor = null;
    ui.empty(host);
    host.removeAttribute('data-mode');
    // The verdict stays cached (the text hasn't changed), but nothing is
    // showing its message anymore — don't mark an editor that isn't there.
    host.removeAttribute('data-invalid');
  };

  /** The fallback: a normal string input plus one line saying what is missing.
   *  Deliberately still editable — a flow loaded from disk carries a condition
   *  that must remain visible and fixable without running anything. */
  const renderPlain = (notice: string, action?: {label: string; run: () => void}): void => {
    reset();
    host.setAttribute('data-mode', 'plain');
    const input = ui.input.string(label, {
      value,
      onValueChanged: (v) => emit(String(v ?? '')),
    });
    const hint = ui.divText(notice, 'ff-expression-editor-hint');
    const rows: HTMLElement[] = [input.root, hint];
    if (action) {
      const button = ui.link(action.label, action.run, 'Run the flow up to the upstream node');
      button.classList.add('ff-expression-editor-load');
      rows.push(button);
    }
    host.appendChild(ui.divV(rows));
  };

  const renderFormula = async (table: DG.DataFrame): Promise<void> => {
    const func = DG.Func.find(EDITOR_FUNC)[0];
    if (!func)
      return renderPlain('Install PowerPack to edit this with the formula editor.');

    // The editor is driven by an AddNewColumn call: it reads `table` for
    // autocomplete and `type` for validation, and publishes the edited text
    // back onto `expression`. Nothing is ever added to the table — both
    // `expressionEditorOnly` and (for conditions) `filterFormulaEditor` make
    // the platform's own AddNewColumn return an empty column list.
    const call = DG.Func.byName('AddNewColumn').prepare({
      table, expression: value, name: label,
      type: options.booleanOnly ? DG.COLUMN_TYPE.BOOL : 'auto',
    });
    call.setAuxValue('expressionEditorOnly', true);
    if (options.booleanOnly) call.setAuxValue('filterFormulaEditor', true);

    let widget: DG.Widget;
    try {
      widget = await func.apply({call}) as DG.Widget;
    } catch (e) {
      console.error('Flow: could not open the formula editor', e);
      return renderPlain('The formula editor could not be opened — edit the expression as text.');
    }
    // A newer render may have replaced us while the widget was being built.
    if (builtFor !== table) return;

    // Read the value back OFF THE CALL, not off the event payload:
    // `FuncCallParam.onChanged` emits the changed parameter, not its value, so
    // stringifying the payload stores a literal "[object Object]" — which then
    // reaches the function as the condition and fails at run time.
    const subscription = call.inputParams['expression'].onChanged.subscribe(() => {
      let next = '';
      try {
        next = String(call.getParamValue('expression') ?? '');
      } catch {/* the call is going away */}
      emit(next);
    });

    disposeWidget = (): void => {
      subscription.unsubscribe();
      try {
        widget.detach?.();
      } catch {/* the widget is going away either way */}
    };

    ui.empty(host);
    host.setAttribute('data-mode', 'formula');
    // Label BESIDE the editor, not above it: the plain fallback is a real DG
    // input and every neighbouring row uses the panel's 90px label column, so a
    // stacked label made this the one row that broke the form grid. (The
    // platform does the same for its own multi-line inputs — see the panel's
    // Description textarea.)
    host.appendChild(ui.divH([
      ui.divText(label, 'ff-expression-editor-label'),
      widget.root,
    ]));
  };

  const build = (): void => {
    const table = ctx.table?.(options.tableParam) ?? null;
    if (table) {
      if (builtFor === table) return; // already showing this exact table
      reset();
      builtFor = table;
      void renderFormula(table);
      return;
    }
    if (!ctx.isConnected?.(options.tableParam)) {
      renderPlain(`Connect a table to “${options.tableParam}” to pick columns and validate as you type.`);
      return;
    }
    renderPlain('The upstream table has not been computed yet — run it to pick columns and validate as you type.', {
      label: 'Run to load',
      run: () => {
        if (!ctx.produceTable) return;
        void ctx.produceTable(options.tableParam)
          .then((produced) => {
            if (produced) build();
            else grok.shell.error('The flow ran but produced no table for this input.');
          })
          .catch((e) => grok.shell.error(`Could not load the table: ${e?.message ?? e}`));
      },
    });
  };

  // Listened for on the host, not on the widget's own root: the widget is
  // mounted inside `host` and its verdicts bubble, so the subscription survives
  // the editor being rebuilt for a different table.
  const onValidated = (e: Event): void => recordVerdict((e as CustomEvent).detail);
  host.addEventListener(VALIDATED_EVENT, onValidated);

  ed.element = host;
  ed.getValue = (): unknown => value;
  ed.setValue = (v): void => {
    value = v === undefined || v === null ? '' : String(v);
    build();
  };
  ed.detach = (): void => {
    host.removeEventListener(VALIDATED_EVENT, onValidated);
    reset();
  };
  return ed;
}
