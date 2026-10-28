import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {LoadedModel, loadedModelOf, targetSchemaOf} from '../apply/apply-model';
import {ColumnMapping, compatibility, exactMapping, suggestMapping} from '../apply/column-matching';
import {applicableTables} from '../catalog/applicable-tables';
import {EngineRegistry} from '../engines/engine-registry';
import {ForgeError} from '../forge-error';
import {ModelRow} from '../generated/db';
import {evaluateModel, storedTrainingTable} from '../training/evaluation';
import {ButtonGate} from './button-gate';
import {storedSummary} from './metrics-table';
import {reportError} from './report-error';
import {ResultsPane} from './results-pane';

const NO_TABLE = 'No open table has the columns the model needs. Open its training data or a similar table.';

/** **Run Evaluation**: the model measured on the data it keeps a link to or a copy of, else on an open table the user
 * chooses in **Evaluate model**. The data loaded for it is not added to the workspace. */
export async function runEvaluation(model: ModelRow): Promise<void> {
  try {
    const loaded = loadedModelOf(model, EngineRegistry.discover());
    const table = await storedTrainingTable(model);
    if (table === null) {
      const tables = applicableTables(model, grok.shell.tables);
      if (tables.length === 0)
        throw new ForgeError(NO_TABLE);
      evaluateModelDialog(model, loaded, tables).show();
      return;
    }
    const target = table.col(model.target_name);
    if (target === null) {
      throw new ForgeError(`The data source has no column '${model.target_name}' the model predicts. ` +
        'Choose another table.');
    }
    await EvaluationView.open(model, loaded, table, target, exactMapping(loaded.features, table));
  } catch (e) {
    reportError(e);
  }
}

/** **Evaluate model**: **Table**, one of the open [tables] the model fits, and **Target**, its column with the actual
 * values; OK opens the evaluation. */
export function evaluateModelDialog(model: ModelRow, loaded: LoadedModel, tables: DG.DataFrame[]): DG.Dialog {
  const schema = targetSchemaOf(model.target) ??
    {name: model.target_name, type: model.task === 'regression' ? 'double' : 'string'};
  const filter = (c: DG.Column) => compatibility(schema, c, loaded.engine.name).kind !== 'error';
  const presetOf = (t: DG.DataFrame) => {
    const col = t.col(model.target_name);
    return col !== null && filter(col) ? col : null;
  };
  const current = grok.shell.currentTable;
  const table = tables.find((t) => t.dart === current?.dart) ?? tables[0];
  const targetInput = ui.input.column('Target', {table, filter, value: presetOf(table) ?? undefined, nullable: false,
    tooltipText: 'Column with the actual values the predictions are compared with.',
    onValueChanged: () => gate.update()});
  const tableInput = ui.input.table('Table', {items: tables, value: table, nullable: false,
    tooltipText: 'Table to evaluate the model on; only the open tables with the columns the model needs are listed.',
    onValueChanged: (t) => {
      if (t !== null) {
        ui.input.setColumnInputTable(targetInput, t, filter);
        targetInput.value = presetOf(t);
      }
      gate.update();
    }});
  const dialog = ui.dialog('Evaluate model').add(tableInput).add(targetInput);
  dialog.onOK(() => {
    const [t, target] = [tableInput.value, targetInput.value];
    if (t !== null && target !== null)
      void EvaluationView.open(model, loaded, t, target, suggestMapping(loaded.features, t));
  });
  const gate = new ButtonGate(dialog.getButton('OK'),
    () => tableInput.value === null ? 'Choose a table.' : targetInput.value === null ? 'Choose a target.' : null);
  gate.update();
  return dialog;
}

/** **Evaluation: <model>**: the model's results on the evaluated table in a {@link ResultsPane}, the stored metrics
 * with the **Evaluation** column next to them. The view holds the evaluated table; nothing is added to it. */
export class EvaluationView extends DG.ViewBase {
  readonly resultsPane = new ResultsPane();
  /** The results; `aria-busy` until they are shown, for automation. */
  private readonly results = ui.panel([this.resultsPane.root], 'forge-train-results');
  private isClosed = false;

  private constructor(model: ModelRow) {
    super();
    this.name = `Evaluation: ${model.name}`;
    this.box = true;
    this.results.setAttribute('aria-busy', 'true');
    this.root.appendChild(ui.box(this.results));
  }

  getIcon(): HTMLElement {
    return ui.iconSvg('model');
  }

  detach(): void {
    this.isClosed = true;
    this.resultsPane.clear();
    super.detach();
  }

  /** Adds the view and evaluates [model] on [table] under the cancellable task-bar progress `Evaluating <model>`; a
   * cancel or a failure closes the view. */
  static async open(model: ModelRow, loaded: LoadedModel, table: DG.DataFrame, target: DG.Column,
    mapping: ColumnMapping): Promise<EvaluationView> {
    const view = new EvaluationView(model);
    grok.shell.addView(view);
    view.resultsPane.showHint('Evaluating the model...');
    const stored = storedSummary(model) ??
      {metrics: {train: {}, validation: {}}, rowCount: model.row_count ?? 0, skippedRows: 0};
    const progress = DG.TaskBarProgressIndicator.create(`Evaluating ${model.name}`, {cancelable: true});
    try {
      const result = await evaluateModel({model: loaded, table, target, mapping,
        positiveClass: stored.metrics.positiveClass}, progress);
      if (view.isClosed)
        return view;
      const {values, rowCount, skippedRows, positiveClass} = result;
      await view.resultsPane.show({engine: loaded.engine, table, columns: result.columns, target: result.target,
        prediction: result.prediction, score: result.score, curves: result.curves,
        keptRows: result.keptRows, positiveClass,
        cutoff: loaded.options.binaryClassificationThreshold, blob: result.blob,
        summary: {...stored, evaluation: {values, tableName: table.name, rowCount, skippedRows,
          ...(result.classes === undefined ? {} : {classes: result.classes})}}});
    } catch (e) {
      if (progress.canceled)
        grok.shell.warning('Evaluation was cancelled.');
      else
        reportError(e);
      view.close();
    } finally {
      progress.close();
      view.results.setAttribute('aria-busy', 'false');
    }
    return view;
  }
}
