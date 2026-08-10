/** Formula editor for expression-holding string parameters — mounts PowerPack's Add New Column
 *  editor inline; falls back to a plain string input when no table (or no PowerPack) is available. */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import type {
  CustomEditorContext, CustomInputEditor, CustomInputEditorFactory, FuncNodeValidator,
} from '../../utils/func-input-overrides';
import {getParamDisplayName} from '../../utils/dart-proxy-utils';
import {propertyNameToFriendly} from '../../utils/naming';

/** Absent PowerPack → the plain string input, which is exactly the no-table state. */
const EDITOR_FUNC = {package: 'PowerPack', name: 'expressionEditorWidget'};

/** PowerPack's validity contract, duplicated as a literal — packages don't import each other's source. */
const VALIDATED_EVENT = 'expression-validated';

/** Last verdict per node — weakly held, deliberately NOT a node property (serializing it would dirty a flow nobody edited). */
const validated = new WeakMap<object, {expression: string; error: string}>();

/** Readiness check for an expression param. Fails OPEN when it cannot decide — a flow loaded from disk and never opened must not be permanently unrunnable. */
export function expressionRequirements(paramName: string): FuncNodeValidator {
  const label = propertyNameToFriendly(paramName);
  return (node) => {
    const verdict = validated.get(node);
    if (!verdict || verdict.error === '') return [];
    if (verdict.expression !== String(node.inputValues[paramName] ?? '')) return [];
    // First sentence only — it has to fit on a node's hint line.
    return [`${label} — ${verdict.error.split('. ')[0]}`];
  };
}

export interface ExpressionEditorOptions {
  tableParam: string;
  /** Constrain to a true/false formula (`aux.filterFormulaEditor`). */
  booleanOnly: boolean;
}

export function expressionEditor(options: ExpressionEditorOptions): CustomInputEditorFactory {
  return (param, ctx) => buildExpressionEditor(param, ctx, options);
}

export const rowConditionEditor = expressionEditor({tableParam: 'table', booleanOnly: true});

export const columnFormulaEditor = expressionEditor({tableParam: 'table', booleanOnly: false});

function buildExpressionEditor(
  param: DG.Property, ctx: CustomEditorContext, options: ExpressionEditorOptions,
): CustomInputEditor {
  const ed: CustomInputEditor = {} as CustomInputEditor;
  const host = ui.div([], 'ff-expression-editor');
  const label = getParamDisplayName(param);
  let value = '';
  /** The table the mounted editor was built against — rebuilding on every render would throw away the cursor and undo history. */
  let builtFor: DG.DataFrame | null = null;
  let disposeWidget: (() => void) | null = null;

  const emit = (next: string): void => {
    if (next === value) return;
    value = next;
    ed.onChanged?.(next);
  };

  const recordVerdict = (detail: {expression: string; error: string}): void => {
    const error = String(detail?.error ?? '');
    host.toggleAttribute('data-invalid', error !== '');
    const node = ctx.node;
    if (!node) return;
    // "No verdict yet" counts as CLEAN: every mount validates, so treating the first clean answer as a change would report a parameter edit on merely opening the panel.
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
    // The verdict stays cached, but nothing shows its message anymore — don't mark an editor that isn't there.
    host.removeAttribute('data-invalid');
  };

  /** Fallback string input — deliberately still editable, so a flow loaded from disk stays fixable without running anything. */
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

    // Both `expressionEditorOnly` and (for conditions) `filterFormulaEditor` make the platform's AddNewColumn add nothing to the table.
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

    // Read the value back OFF THE CALL — `FuncCallParam.onChanged` emits the parameter, not its value; stringifying the payload stores "[object Object]".
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
    // Label BESIDE the editor — a stacked label broke the panel's 90px-label form grid.
    host.appendChild(ui.divH([
      ui.divText(label, 'ff-expression-editor-label'),
      widget.root,
    ]));
  };

  const build = (): void => {
    const table = ctx.table?.(options.tableParam) ?? null;
    if (table) {
      if (builtFor === table) return;
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

  // Listen on the host, not the widget root — verdicts bubble, so the subscription survives editor rebuilds.
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
