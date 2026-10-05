import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {randomInt} from '@datagrok-libraries/utils/src/random';
import {defaultHyperparameters, Engine, Hyperparameters, hyperparametersOf, isComplete} from '../engines/engine';
import {EngineRegistry} from '../engines/engine-registry';
import {errorMessage, ForgeError} from '../forge-error';
import {TrainingRunStatus} from '../generated/db';
import {METRIC_DESCRIPTIONS, METRIC_IDS, METRIC_LABELS, MetricId} from '../metrics/metrics';
import {DatasetFingerprint, datasetFingerprint} from '../storage/dataset-fingerprint';
import {modelFieldsOf, trainingRunOf} from '../storage/model-fields';
import {saveModel} from '../storage/model-store';
import {linkTrainingRun, recordTrainingRun} from '../storage/training-run-store';
import {missingColumnsOf} from '../preparation/missing-values';
import {releaseFrame} from '../preparation/shared-frame';
import {defaultFeatures} from '../training/default-features';
import {checkTrainable, MetricsRecord, prepareTraining, TrainingProblems, trainingProblems, TrainingRequest,
  TrainingResult, TrainingSelection, trainModel} from '../training/train-model';
import {ButtonGate} from './button-gate';
import {CollapsibleGroup} from './collapsible-group';
import {MissingValuesInputs} from './missing-values-inputs';
import {reportError} from './report-error';
import {saveModelDialog} from './save-model-dialog';

const TRAIN_ENGINE = 'XGBoost';
const FOLDS = 5;
const CHECK_DELAY_MS = 300;
const NAME_FEATURES_LENGTH = 30;
const RESULTS_HINT = 'Choose the target and the features, then click Train.';
const CHECKING = 'Checking the selection...';
const TRAINING = 'Training is in progress.';
const NO_TARGET = 'Choose a target.';
const TARGET_TOOLTIP = 'The column the model learns to predict.';

interface TrainForm {
  table: DG.DataFrame;
  tableInput: DG.InputBase<DG.DataFrame | null>;
  target: DG.InputBase<DG.Column | null>;
  features: DG.InputBase<DG.Column[]>;
  missingValues: MissingValuesInputs;
  hyperparameters: Map<string, DG.InputBase>;
  groups: CollapsibleGroup[];
}

export interface Training {
  result: TrainingResult;
  runId: string;
  datasetName: string;
  fingerprint: DatasetFingerprint;
  isSaved: boolean;
}

export class TrainView extends DG.ViewBase {
  readonly trainButton: HTMLButtonElement;
  readonly saveButton: HTMLButtonElement;
  lastTraining: Training | undefined;
  private readonly trainGate: ButtonGate;
  private readonly engine: Engine;
  private readonly formPane: HTMLDivElement;
  private readonly resultsHost: HTMLDivElement;
  private readonly changes = new rxjs.Subject<void>();
  private form: TrainForm;
  // The current form's subscriptions: a Table change rebuilds the form and drops the old ones.
  private formSubs: rxjs.Subscription[] = [];
  private problems: TrainingProblems | undefined;
  private checkNumber = 0;
  private trainBlocker: string | null = null;
  private isTraining = false;

  private constructor(table: DG.DataFrame, engine: Engine) {
    super();
    this.name = 'Predictive model';
    this.box = true;
    this.engine = engine;
    this.trainButton = ui.bigButton('Train', () => this.train());
    this.trainGate = new ButtonGate(this.trainButton, () => this.isTraining ? TRAINING : this.trainBlocker);
    this.saveButton = ui.bigButton('Save', () => this.save());
    this.setRibbonPanels([[this.saveButton]]);
    this.formPane = ui.panel([]);
    this.resultsHost = ui.div([]);
    this.subs.push(DG.debounce(this.changes, CHECK_DELAY_MS).subscribe(() => this.revalidate()));
    this.form = this.createForm(table);
    this.root.appendChild(ui.splitH([this.formPane, ui.panel([ui.h2('Results'), this.resultsHost])]));
  }

  getIcon(): HTMLElement {
    return ui.iconSvg('model');
  }

  detach(): void {
    this.unsubscribeForm();
    super.detach();
  }

  get tableInput(): DG.InputBase<DG.DataFrame | null> {
    return this.form.tableInput;
  }

  get targetInput(): DG.InputBase<DG.Column | null> {
    return this.form.target;
  }

  get featuresInput(): DG.InputBase<DG.Column[]> {
    return this.form.features;
  }

  get hyperparameterInputs(): ReadonlyMap<string, DG.InputBase> {
    return this.form.hyperparameters;
  }

  get missingValuesInputs(): MissingValuesInputs {
    return this.form.missingValues;
  }

  /** **Data** and **Method**. */
  get groups(): readonly CollapsibleGroup[] {
    return this.form.groups;
  }

  static async create(table: DG.DataFrame): Promise<TrainView> {
    const engine = EngineRegistry.discover().find((e) => e.name === TRAIN_ENGINE && isComplete(e));
    if (engine === undefined)
      throw new ForgeError(`The ${TRAIN_ENGINE} method is not installed. Install the EDA package.`);
    return new TrainView(table, engine);
  }

  static async open(): Promise<void> {
    try {
      const table = grok.shell.currentTable;
      if (table === null)
        grok.shell.warning('Open a table first.');
      else
        grok.shell.addView(await TrainView.create(table));
    } catch (e) {
      reportError(e);
    }
  }

  async train(): Promise<void> {
    if (this.isTraining)
      return;
    this.isTraining = true;
    this.updateButtons();
    let progress: DG.TaskBarProgressIndicator | undefined;
    try {
      const selection = this.selection(randomInt(2 ** 31));
      if (selection === null)
        throw new ForgeError(NO_TARGET);
      await checkTrainable(selection);
      this.clearResults();
      progress = DG.TaskBarProgressIndicator.create(`Training ${selection.engine.name} model`, {cancelable: true});
      const request = await prepareTraining(selection);
      try {
        this.lastTraining = await this.trainAndRecord(request, progress);
      } finally {
        releaseFrame(request.features);
      }
      this.showResults(this.lastTraining);
    } catch (e) {
      reportError(e);
    } finally {
      progress?.close();
      this.isTraining = false;
      this.updateButtons();
    }
  }

  save(): void {
    const training = this.lastTraining;
    if (training === undefined || training.isSaved)
      return;
    saveModelDialog(TrainView.defaultModelName(training.result),
      (name, description) => this.saveModelAs(name, description)).show();
  }

  async saveModelAs(name: string, description: string): Promise<void> {
    const training = this.lastTraining;
    if (training === undefined || training.isSaved)
      return;
    training.isSaved = true;
    this.updateButtons();
    const {result, runId, datasetName, fingerprint} = training;
    let id: string;
    try {
      id = await saveModel(modelFieldsOf({name, description, engine: this.engine, datasetName, result, fingerprint}),
        result.blob);
    } catch (e) {
      training.isSaved = false;
      this.updateButtons();
      throw e;
    }
    await linkTrainingRun(runId, id);
    grok.shell.info(`Model "${name}" saved. Open ML | Forge | Models to see it.`);
  }

  private createForm(table: DG.DataFrame): TrainForm {
    const target = table.columns.byIndex(table.columns.length - 1);
    const onChanged = () => this.requestCheck();
    const targetInput = ui.input.column('Target', {table, value: target, nullable: false,
      tooltipText: TrainView.targetTooltip(target), onValueChanged: (t, input) => {
        input.setTooltip(TrainView.targetTooltip(t));
        onChanged();
      }});
    const checked = defaultFeatures(table, target).map((c) => c.name);
    const missingValues = new MissingValuesInputs(onChanged,
      () => TrainView.messageOf(this.problems?.missingValues));
    const featuresInput = ui.input.columns('Features', {table, checked, nullable: false,
      tooltipText: 'Columns the model uses to make predictions.', onValueChanged: (columns) => {
        missingValues.update(missingColumnsOf(columns));
        onChanged();
      }});
    missingValues.update(missingColumnsOf(featuresInput.value));
    targetInput.addValidator(() => TrainView.messageOf(this.problems?.target));
    featuresInput.addValidator(() => TrainView.messageOf(this.problems?.features));
    const tableInput = ui.input.table('Table', {items: grok.shell.tables, value: table,
      tooltipText: 'Data to train the model on.', onValueChanged: (t) => {
        if (t !== null)
          this.form = this.createForm(t);
      }});
    const engineName = this.engine.name;
    const methodInput = ui.input.choice('Method', {items: [engineName], value: engineName, nullable: false,
      tooltipText: 'The machine learning method that builds the model.'});

    const defaults = defaultHyperparameters(this.engine);
    const hyperparameters = new Map<string, DG.InputBase>();
    for (const p of hyperparametersOf(this.engine)) {
      const input = ui.input.forProperty(p);
      const value = defaults[p.name];
      if (value !== undefined)
        input.value = value;
      hyperparameters.set(p.name, input);
    }

    const dataInputs = [tableInput, targetInput, featuresInput, ...missingValues.inputs];
    const methodInputs = [methodInput, ...hyperparameters.values()];
    const groups = [new CollapsibleGroup('Data', [ui.form(dataInputs)]),
      new CollapsibleGroup('Method', [ui.form(methodInputs)])];
    this.unsubscribeForm();
    this.formSubs = [...groups[0].expandOnError(dataInputs), ...groups[1].expandOnError(methodInputs)];
    const trainRow = ui.buttonsInput([this.trainButton]);
    trainRow.classList.add('forge-train-row');
    const buttons = ui.form([]);
    buttons.append(trainRow);
    ui.empty(this.formPane);
    this.formPane.append(...groups.map((g) => g.root), buttons);
    this.clearResults();
    this.requestCheck();
    return {table, tableInput, target: targetInput, features: featuresInput, missingValues, hyperparameters, groups};
  }

  private unsubscribeForm(): void {
    for (const sub of this.formSubs)
      sub.unsubscribe();
    this.formSubs = [];
  }

  private selection(seed: number): TrainingSelection | null {
    const target = this.form.target.value;
    if (target === null)
      return null;
    const hyperparameters: Hyperparameters = {};
    for (const [name, input] of this.form.hyperparameters) {
      const value: unknown = input.value;
      if (typeof value === 'number' || typeof value === 'string' || typeof value === 'boolean')
        hyperparameters[name] = value;
    }
    return {engine: this.engine, features: this.form.features.value, target, hyperparameters, seed, folds: FOLDS,
      missingValues: this.form.missingValues.settings()};
  }

  private requestCheck(): void {
    this.checkNumber++;
    this.problems = undefined;
    this.trainBlocker = CHECKING;
    ui.setDisabled(this.trainButton, true);
    this.changes.next();
  }

  private async revalidate(): Promise<void> {
    const checkNumber = this.checkNumber;
    try {
      const selection = this.selection(0);
      const problems = selection === null ? {target: [NO_TARGET], features: [], missingValues: []} :
        await trainingProblems(selection);
      if (checkNumber !== this.checkNumber)
        return;
      this.problems = problems;
      const settings = this.form.missingValues.visibleInputs;
      for (const input of [this.form.target, this.form.features, ...settings])
        input.validate();
      const invalidSettings = settings.map((input) => input.validity).filter((v): v is string => v !== null);
      this.trainBlocker = [...problems.target, ...problems.features, ...problems.missingValues,
        ...invalidSettings][0] ?? null;
      this.updateButtons();
    } catch (e) {
      if (checkNumber !== this.checkNumber)
        return;
      this.trainBlocker = errorMessage(e);
      this.updateButtons();
      reportError(e, true);
    }
  }

  private updateButtons(): void {
    this.trainGate.update();
    const training = this.lastTraining;
    ui.setDisabled(this.saveButton, this.isTraining || training === undefined || training.isSaved);
  }

  private async trainAndRecord(request: TrainingRequest, progress: DG.TaskBarProgressIndicator): Promise<Training> {
    const datasetName = this.form.table.name;
    const fingerprint = datasetFingerprint(request.features, request.target);
    const startedOn = new Date();
    const record = (status: TrainingRunStatus, metrics?: MetricsRecord, error?: string) =>
      recordTrainingRun(trainingRunOf({request, datasetName, fingerprint, status, startedOn: startedOn.toISOString(),
        durationMs: Date.now() - startedOn.getTime(), metrics, error}));
    let result: TrainingResult;
    try {
      result = await trainModel(request, progress);
    } catch (e) {
      if (progress.canceled)
        await record('cancelled');
      else
        await record('failed', undefined, errorMessage(e));
      throw e;
    }
    const runId = await record('completed', result.metrics);
    return {result, runId, datasetName, fingerprint, isSaved: false};
  }

  private clearResults(): void {
    this.lastTraining = undefined;
    ui.setDisabled(this.saveButton, true);
    ui.empty(this.resultsHost);
    this.resultsHost.append(ui.divText(RESULTS_HINT));
  }

  private showResults({result}: Training): void {
    const {train, validation, positiveClass} = result.metrics;
    const ids = METRIC_IDS.filter((id) => validation[id] !== undefined);
    const label = (id: MetricId) => ui.tooltip.bind(ui.label(METRIC_LABELS[id]), METRIC_DESCRIPTIONS[id]);
    const format = (value: number | undefined) => value === undefined ? '' : value.toFixed(3);
    const skippedRows = result.options.missingValues?.skippedRows ?? 0;
    const rows = skippedRows === 0 ? [] :
      [ui.divText(`Rows: ${result.rowCount} used, ${skippedRows} skipped (missing values).`)];
    ui.empty(this.resultsHost);
    this.resultsHost.append(
      ui.table(ids, (id) => [label(id), format(train[id]), format(validation[id])], ['Metric', 'Train', 'Validation']),
      ...rows,
      ui.divText(`Validation: ${FOLDS}-fold cross-validation on ${result.rowCount} rows, seed ${result.seed}.`),
      ...(positiveClass === undefined ? [] : [ui.divText(`Positive class: ${positiveClass}.`)]),
    );
  }

  private static defaultModelName(result: TrainingResult): string {
    const features = result.features.columns.map((c) => c.name).join(', ');
    const shortened = features.length > NAME_FEATURES_LENGTH ? `${features.substring(0, NAME_FEATURES_LENGTH)}...` :
      features;
    return `Predict ${result.target.name} by ${shortened}`;
  }

  private static targetTooltip(target: DG.Column | null): string {
    const count = target?.stats.missingValueCount ?? 0;
    if (target === null || count === 0)
      return TARGET_TOOLTIP;
    const skipped = count === 1 ? '1 row without a value is skipped' : `${count} rows without a value are skipped`;
    return `${TARGET_TOOLTIP} ${target.name}: ${skipped}.`;
  }

  private static messageOf(messages: string[] | undefined): string | null {
    return messages?.length ? messages.join(' ') : null;
  }
}
