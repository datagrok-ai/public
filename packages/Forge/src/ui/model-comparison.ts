import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {compareFormFields, compareModels, CompareModelRow} from '../catalog/compare-models';
import {isRecord} from '../preparation/preparation-options';
import {readOnlyGrid} from './data-grid';
import {modelIcon} from './model-panes';

const COMPARISON_TYPE = 'forge.model.comparison';
const HANDLER_NAME = 'Forge model comparison handler';
const FORMS = 'Forms';
const CARD_HEIGHT = 160;
const MAX_SHOWN_CARDS = 4;
export const NO_FORMS = 'Install the PowerGrid package to see the comparison as forms.';

/** Two or more models chosen together: the context panel shows them side by side. */
export class ModelComparison {
  readonly kind = COMPARISON_TYPE;

  constructor(readonly rows: CompareModelRow[]) {}
}

/** Shows a {@link ModelComparison} in the context panel as forms, one card per model. */
export class ModelComparisonHandler extends DG.ObjectHandler<ModelComparison> {
  get type(): string {
    return COMPARISON_TYPE;
  }

  get name(): string {
    return HANDLER_NAME;
  }

  /** Registers the handler once per page, unless [registered] (the handlers' names) has it already. */
  static registerOnce(registered: string[]): void {
    if (!registered.includes(HANDLER_NAME))
      DG.ObjectHandler.register(new ModelComparisonHandler());
  }

  // By its kind, not its class: the package-test bundle builds comparisons of its own copy of the class.
  isApplicable(x: unknown): boolean {
    return isRecord(x) && x.kind === COMPARISON_TYPE && Array.isArray(x.rows);
  }

  getCaption(x: ModelComparison): string {
    return captionOf(x.rows.length);
  }

  renderProperties(x: ModelComparison): HTMLElement {
    return comparisonForms(x.rows);
  }
}

export function hasFormsViewer(): boolean {
  return DG.Func.find({package: 'PowerGrid', name: 'formsViewer'}).length > 0;
}

/** The Forms viewer over [df]: a card per row, with the summary and the validation metrics. */
function comparisonFormsViewer(df: DG.DataFrame): DG.Viewer {
  df.selection.setAll(true);
  return DG.Viewer.fromType(FORMS, df, {fieldsColumnNames: compareFormFields(df), showCurrentRow: false,
    showMouseOverRow: false, showSelectedRows: true});
}

/** The comparison of [rows] for the context panel: the title `Compare N models` as a model's panel has it, then the
 * Forms viewer, or without PowerGrid the grid and why. */
export function comparisonForms(rows: CompareModelRow[]): HTMLElement {
  const df = compareModels(rows);
  // The context panel shows no caption of its own: the title is the platform's accordion title, as in modelAccordion.
  const title = ui.div([ui.span([modelIcon(), ui.label(captionOf(rows.length))])], 'd4-accordion-title');
  if (!hasFormsViewer())
    return ui.divV([title, ui.divText(NO_FORMS), readOnlyGrid(df).root]);
  const forms = comparisonFormsViewer(df);
  // Inline, computed: the viewer lays its cards out in the size it is given when attached.
  forms.root.style.width = '100%';
  forms.root.style.height = `${Math.min(rows.length, MAX_SHOWN_CARDS) * CARD_HEIGHT}px`;
  return ui.divV([title, forms.root]);
}

function captionOf(count: number): string {
  return `Compare ${count} models`;
}

/** The table view **Compare models** with the grid and, with PowerGrid, the Forms viewer. */
export function openComparisonView(rows: CompareModelRow[]): void {
  const df = compareModels(rows);
  const view = grok.shell.addTableView(df);
  if (hasFormsViewer())
    view.addViewer(comparisonFormsViewer(df));
  else
    grok.shell.info(NO_FORMS);
}
