import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {randomInt} from '@datagrok-libraries/utils/src/random';
import {ModelInfo} from '../catalog/model-edit';
import {defaultHyperparameters, Engine, Hyperparameters, hyperparametersOf, isComplete} from '../engines/engine';
import {EngineRegistry} from '../engines/engine-registry';
import {errorMessage, ForgeError} from '../forge-error';
import {TrainingRunStatus} from '../generated/db';
import {_package} from '../package';
import {DatasetFingerprint, datasetFingerprint} from '../storage/dataset-fingerprint';
import {deleteTrainingCopy, uploadTrainingCopy} from '../storage/dataset-copy';
import {datasetRefOf} from '../storage/dataset-ref';
import {ModelStorage, modelFieldsOf, TrainingRunRecord, trainingRunOf} from '../storage/model-fields';
import {saveModel} from '../storage/model-store';
import {linkTrainingRun, recordTrainingRun} from '../storage/training-run-store';
import {missingColumnsOf} from '../preparation/missing-values';
import {releaseFrame} from '../preparation/shared-frame';
import {defaultFeatures} from '../training/default-features';
import {checkSelection, failedCheck, hasDataProblems, MetricsRecord, prepareTraining, retrainsLive, SelectionCheck,
  TrainingResult, TrainingSelection, trainModel} from '../training/train-model';
import {TrainingQueue} from '../training/training-queue';
import {openApplyDialog} from './apply-model-dialog';
import {ButtonGate} from './button-gate';
import {CollapsibleGroup} from './collapsible-group';
import {ForgeApp} from './forge-app';
import {MissingValuesInputs} from './missing-values-inputs';
import {metricsTable} from './model-panes';
import {reportError} from './report-error';
import {saveModelDialog, StorageChoice} from './save-model-dialog';

const PREFERRED_ENGINE = 'XGBoost';
const FOLDS = 5;
/** The built-in tool's delay after the last change before the preview is trained again. */
export const CHECK_DELAY_MS = 200;
const NAME_FEATURES_LENGTH = 30;
const RESULTS_HINT = 'Choose the target and the features; the model trains as you change them.';
const FIX_SETTINGS = 'Fix the settings.';
const FORM_WIDTH = 400;
const TRAIN_TOOLTIP = 'Train the model on the current selection.';
const CHECKING = 'Checking the selection...';
const TRAINING = 'Training is in progress.';
const TRAINED = 'The model is trained on this selection. Change an input to train again.';
const CANCELLED = 'Training was cancelled.';
const NO_TARGET = 'Choose a target.';
const TARGET_TOOLTIP = 'The column the model learns to predict.';

interface TrainForm {
  table: DG.DataFrame;
  tableInput: DG.InputBase<DG.DataFrame | null>;
  target: DG.InputBase<DG.Column | null>;
  features: DG.InputBase<DG.Column[]>;
  missingValues: MissingValuesInputs;
  method: DG.ChoiceInput<string | null>;
  hyperparameters: Map<string, DG.InputBase>;
  /** The Method group's content: rebuilt with the method. */
  methodBody: HTMLDivElement;
  groups: CollapsibleGroup[];
}

export interface Training {
  result: TrainingResult;
  runId: string;
  engine: Engine;
  table: DG.DataFrame;
  /** The user's feature and target columns the model was trained on: what a copy uploads. */
  columns: DG.Column[];
  datasetName: string;
  fingerprint: DatasetFingerprint;
  isSaved: boolean;
}

export class TrainView extends DG.ViewBase {
  /** **Train** under the inputs: shown only for a method that does not retrain on every change. */
  readonly trainButton: HTMLButtonElement;
  readonly saveButton: HTMLButtonElement;
  lastTraining: Training | undefined;
  private readonly trainRow: HTMLDivElement;
  private readonly trainGate: ButtonGate;
  /** The complete methods discovered when the view opened, in discovery order. */
  private readonly engines: Engine[];
  private readonly formPane: HTMLDivElement;
  private readonly resultsStatus = ui.div([], 'forge-results-status');
  private readonly resultsBody = ui.div([]);
  private readonly changes = new rxjs.Subject<void>();
  private readonly queue = new TrainingQueue();
  /** Hyperparameter values per method name, kept for the session of the view (a Table change keeps them). */
  private readonly hyperparameterValues = new Map<string, Hyperparameters>();
  /** Methods whose check failed, logged once per view. */
  private readonly loggedFailures = new Set<string>();
  /** Whether each method retrains on every change, as the last check of this table that could ask said. */
  private readonly interactivity = new Map<string, boolean>();
  private engine: Engine;
  /** The method the user picked; null while the suggested one applies. */
  private userChoice: string | null = null;
  private form: TrainForm;
  // The current form's subscriptions: a Table change rebuilds the form and drops the old ones.
  private formSubs: rxjs.Subscription[] = [];
  // The hyperparameter inputs' subscriptions: a Method change rebuilds them.
  private methodSubs: rxjs.Subscription[] = [];
  /** The latest check's result; undefined while it is pending. */
  private check: SelectionCheck | undefined;
  private checkNumber = 0;
  private activeTrainings = 0;
  /** A training of the current selection completed or failed; any change clears it. */
  private isTrained = false;
  private isSettingMethod = false;

  private constructor(table: DG.DataFrame, engines: Engine[], engine: Engine) {
    super();
    this.name = 'Predictive model';
    this.box = true;
    this.engines = engines;
    this.engine = engine;
    this.trainButton = ui.button('Train', () => void this.train());
    this.trainRow = ui.buttonsInput([this.trainButton]);
    this.trainRow.classList.add('forge-train-row');
    this.trainGate = new ButtonGate(this.trainButton, () => this.blocker(), TRAIN_TOOLTIP);
    this.saveButton = ui.bigButton('Save', () => this.save());
    this.setRibbonPanels([[this.saveButton]]);
    this.formPane = ui.panel([], 'forge-train-form');
    this.subs.push(DG.debounce(this.changes, CHECK_DELAY_MS).subscribe(() => this.revalidate()));
    this.form = this.createForm(table);
    const results = ui.panel([ui.divH([ui.h2('Results'), this.resultsStatus], 'forge-results-header'),
      this.resultsBody], 'forge-train-results');
    const left = ui.box(this.formPane);
    left.style.width = `${FORM_WIDTH}px`;
    this.root.appendChild(ui.splitH([left, ui.box(results)], {}, true));
  }

  getIcon(): HTMLElement {
    return ui.iconSvg('model');
  }

  detach(): void {
    this.queue.supersede();
    this.unsubscribeForm();
    this.unsubscribeMethod();
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

  get methodInput(): DG.ChoiceInput<string | null> {
    return this.form.method;
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

  /** A training is queued, running or being recorded. */
  get isTraining(): boolean {
    return this.activeTrainings > 0;
  }

  static async create(table: DG.DataFrame): Promise<TrainView> {
    const engines = EngineRegistry.discover().filter(isComplete);
    const engine = engines.find((e) => e.name === PREFERRED_ENGINE) ?? engines[0];
    if (engine === undefined)
      throw new ForgeError('No machine learning method is installed. Install the EDA package.');
    return new TrainView(table, engines, engine);
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

  /** **Train**: trains the current selection, whether the method retrains live or not. */
  async train(): Promise<void> {
    if (this.problem() === null && !this.isTraining)
      await this.startTraining(false);
  }

  save(): void {
    const training = this.lastTraining;
    if (training === undefined || training.isSaved)
      return;
    saveModelDialog(TrainView.defaultModelName(training.result), (info, choice) => this.saveModelAs(info, choice),
      {ref: datasetRefOf(training.table), rowCount: training.table.rowCount}).show();
  }

  /** Saves the last training once; `copy` uploads the training columns first, `reference` links the table's origin. */
  async saveModelAs({name, description, tags}: ModelInfo, choice: StorageChoice = {mode: 'none'}): Promise<void> {
    const training = this.lastTraining;
    if (training === undefined || training.isSaved)
      return;
    training.isSaved = true;
    this.updateControls();
    const {result, runId, datasetName, fingerprint, engine, table} = training;
    let tableId: string | null = null;
    let id: string;
    try {
      let storage: ModelStorage = choice.mode === 'reference' ? choice : {mode: 'none'};
      if (choice.mode === 'copy') {
        tableId = await uploadTrainingCopy(training.columns, name);
        storage = {mode: 'copy', tableId};
      }
      id = await saveModel(modelFieldsOf({name, description, tags, engine, datasetName, result, fingerprint, storage}),
        result.blob);
    } catch (e) {
      training.isSaved = false;
      this.updateControls();
      if (tableId !== null)
        await deleteTrainingCopy(tableId).catch((cleanup) => _package.logger.error(errorMessage(cleanup)));
      throw e;
    }
    await linkTrainingRun(runId, id);
    grok.shell.info(ui.divV([ui.divText(`Model "${name}" saved.`), ui.divH([
      ui.link('Apply...', () => void openApplyDialog(table, {modelId: id, preferredTable: table})),
      ui.link('Show in the catalog', () => void ForgeApp.open()),
    ], 'forge-links')]));
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
      () => TrainView.messageOf(this.check?.problems.missingValues));
    const featuresInput = ui.input.columns('Features', {table, checked, nullable: false,
      tooltipText: 'Columns the model uses to make predictions.', onValueChanged: (columns) => {
        missingValues.update(missingColumnsOf(columns));
        onChanged();
      }});
    missingValues.update(missingColumnsOf(featuresInput.value));
    targetInput.addValidator(() => TrainView.messageOf(this.check?.problems.target));
    featuresInput.addValidator(() => TrainView.messageOf(this.check?.problems.features));
    const tableInput = ui.input.table('Table', {items: grok.shell.tables, value: table,
      tooltipText: 'Data to train the model on.', onValueChanged: (t) => {
        if (t !== null)
          this.form = this.createForm(t);
      }});
    const methodInput = ui.input.choice<string | null>('Method', {items: [this.engine.name], value: this.engine.name,
      nullable: false, tooltipText: 'The machine learning method that builds the model.',
      onValueChanged: (name) => this.chooseMethod(name)});
    methodInput.addValidator(() => TrainView.messageOf(this.check?.problems.method));

    const dataInputs = [tableInput, targetInput, featuresInput, ...missingValues.inputs];
    const methodBody = ui.div([]);
    const groups = [new CollapsibleGroup('Data', dataInputs.map((input) => input.root)),
      new CollapsibleGroup('Method', [methodBody])];
    this.unsubscribeForm();
    this.interactivity.clear();
    ui.setDisplay(this.trainRow, false);
    this.formSubs = [...groups[0].expandOnError(dataInputs), ...groups[1].expandOnError([methodInput])];
    ui.empty(this.formPane);
    // One form, so the labels of both groups and the Train row line up.
    const formRoot = ui.form([]);
    formRoot.append(...groups.map((g) => g.root), this.trainRow);
    this.formPane.append(formRoot);
    this.lastTraining = undefined;
    this.showHint(RESULTS_HINT);
    const form: TrainForm = {table, tableInput, target: targetInput, features: featuresInput, missingValues,
      method: methodInput, hyperparameters: new Map(), methodBody, groups};
    this.fillMethod(form);
    this.requestCheck();
    return form;
  }

  /** The Method group of [form] for the current method: Method, its hyperparameters with the values kept for it. */
  private fillMethod(form: TrainForm): void {
    const engine = this.engine;
    const values = {...defaultHyperparameters(engine), ...this.hyperparameterValues.get(engine.name)};
    this.unsubscribeMethod();
    form.hyperparameters.clear();
    for (const p of hyperparametersOf(engine)) {
      const input = ui.input.forProperty(p);
      const value = values[p.name];
      if (value !== undefined)
        input.value = value;
      form.hyperparameters.set(p.name, input);
      this.methodSubs.push(input.onChanged.subscribe(() => {
        this.hyperparameterValues.set(engine.name, TrainView.valuesOf(form.hyperparameters));
        this.requestCheck();
      }));
    }
    const inputs = [...form.hyperparameters.values()];
    this.methodSubs.push(...form.groups[1].expandOnError(inputs));
    ui.empty(form.methodBody);
    form.methodBody.append(form.method.root, ...inputs.map((input) => input.root));
  }

  /** The user picked [name] in **Method**. */
  private chooseMethod(name: string | null): void {
    const engine = this.engines.find((e) => e.name === name);
    if (this.isSettingMethod || engine === undefined)
      return;
    this.userChoice = engine.name;
    if (engine !== this.engine) {
      this.engine = engine;
      this.fillMethod(this.form);
    }
    this.requestCheck();
  }

  private unsubscribeForm(): void {
    for (const sub of this.formSubs)
      sub.unsubscribe();
    this.formSubs = [];
  }

  private unsubscribeMethod(): void {
    for (const sub of this.methodSubs)
      sub.unsubscribe();
    this.methodSubs = [];
  }

  private selection(seed: number): TrainingSelection | null {
    const target = this.form.target.value;
    if (target === null)
      return null;
    return {engine: this.engine, features: this.form.features.value, target,
      hyperparameters: TrainView.valuesOf(this.form.hyperparameters), seed, folds: FOLDS,
      missingValues: this.form.missingValues.settings()};
  }

  /** Every change: the training of the old selection stops, its result goes, and the selection is checked again
   * [CHECK_DELAY_MS] after the last change. */
  private requestCheck(): void {
    this.checkNumber++;
    this.check = undefined;
    this.queue.supersede();
    this.lastTraining = undefined;
    this.isTrained = false;
    this.updateControls();
    this.changes.next();
  }

  private async checkCurrent(): Promise<SelectionCheck> {
    const selection = this.selection(0);
    if (selection === null)
      return failedCheck({target: [NO_TARGET]});
    const check = await checkSelection(selection, this.engines);
    for (const {engine: failed, error} of check.failed) {
      if (!this.loggedFailures.has(failed.name)) {
        this.loggedFailures.add(failed.name);
        _package.logger.error(`The method ${failed.name} could not check the data: ${errorMessage(error)}`);
      }
    }
    return check;
  }

  /** Checks the selection, updates **Method** (the user's choice while it is listed, the suggested method otherwise)
   * and, when the method retrains live, trains. */
  private async revalidate(): Promise<void> {
    const checkNumber = this.checkNumber;
    try {
      let check = await this.checkCurrent();
      if (checkNumber !== this.checkNumber)
        return;
      const engine = this.methodOf(check);
      const target = this.form.target.value;
      if (engine !== this.engine && target !== null) {
        this.engine = engine;
        this.fillMethod(this.form);
        // The data and the list stand; only the new method, which is listed, is asked whether it retrains live.
        const isInteractive = await retrainsLive(engine, this.form.features.value, target);
        if (checkNumber !== this.checkNumber)
          return;
        check = {...check, problems: {...check.problems, method: []}, isInteractive};
      }
      this.check = check;
      this.showMethods(check);
      const settings = this.form.missingValues.visibleInputs;
      const hyperparameters = this.form.hyperparameters.values();
      for (const input of [this.form.target, this.form.features, this.form.method, ...settings, ...hyperparameters])
        input.validate();
      if (check.engines.length > 0)
        this.interactivity.set(this.engine.name, check.isInteractive);
      this.showTrainRow();
      this.updateControls();
      if (this.problem() !== null)
        this.showHint(TrainView.dataProblems(check).length === 0 ? FIX_SETTINGS : RESULTS_HINT);
      else if (check.isInteractive)
        await this.startTraining(true);
      else
        this.showHint(this.notInteractiveHint());
    } catch (e) {
      if (checkNumber !== this.checkNumber)
        return;
      this.check = failedCheck({method: [errorMessage(e)]});
      this.showTrainRow();
      this.updateControls();
      reportError(e, true);
    }
  }

  /** **Train** is shown for a method known not to retrain on every change, never for an interactive or an unknown one:
   * the data checks that fail leave the method's last interactivity on this table. */
  private showTrainRow(): void {
    ui.setDisplay(this.trainRow, this.interactivity.get(this.engine.name) === false);
  }

  /** The method for [check]: the current one while a target or features rule fails or nothing can learn; the user's
   * choice while it is listed; the suggested one otherwise, with a balloon when it replaces the user's choice. */
  private methodOf(check: SelectionCheck): Engine {
    const {problems, engines, best} = check;
    if (hasDataProblems(problems) || engines.length === 0)
      return this.engine;
    const chosen = engines.find((e) => e.name === this.userChoice);
    if (chosen !== undefined)
      return chosen;
    const suggested = best ?? engines[0];
    if (this.userChoice !== null) {
      grok.shell.info(`${this.userChoice} cannot be used with this selection; ${suggested.name} is chosen.`);
      this.userChoice = null;
    }
    return suggested;
  }

  /** **Method** lists the methods that can learn from the selection; it is empty when none can. A failing target or
   * features rule leaves the list as it is. */
  private showMethods({problems, engines}: SelectionCheck): void {
    if (hasDataProblems(problems))
      return;
    this.isSettingMethod = true;
    try {
      this.form.method.items = engines.map((e) => e.name);
      this.form.method.value = engines.includes(this.engine) ? this.engine.name : null;
    } finally {
      this.isSettingMethod = false;
    }
  }

  /** Why **Train** is unavailable: the check is pending, a problem, a training running or the selection already
   * trained; null when it is not. */
  private blocker(): string | null {
    return this.problem() ?? (this.isTraining ? TRAINING : this.isTrained ? TRAINED : null);
  }

  /** The pending check or the selection's first problem, an invalid hyperparameter included; null when the selection
   * can be trained. */
  private problem(): string | null {
    if (this.check === undefined)
      return CHECKING;
    const hyperparameters = [...this.form.hyperparameters.values()].filter((input) => input.validity !== null)
      .map((input) => `${input.caption}: ${input.validity}`);
    const settings = this.form.missingValues.visibleInputs.map((input) => input.validity)
      .filter((v): v is string => v !== null);
    return [...TrainView.dataProblems(this.check), ...settings, ...hyperparameters][0] ?? null;
  }

  /** **Train**, **Save** and the Results hint follow the check, the training and the last result. */
  private updateControls(): void {
    this.trainGate.update();
    const training = this.lastTraining;
    ui.setDisabled(this.saveButton, this.isTraining || training === undefined || training.isSaved);
  }

  /** Trains the current selection through the queue (a newer training supersedes it) with a task-bar progress; records
   * a completed, failed or cancelled run (a superseded one is not recorded) and shows the result. */
  private async startTraining(isLive: boolean): Promise<void> {
    const selection = this.selection(randomInt(2 ** 31));
    if (selection === null)
      return;
    const checkNumber = this.checkNumber;
    const {engine, features, target} = selection;
    const table = this.form.table;
    const datasetName = table.name;
    this.activeTrainings++;
    this.showTraining(isLive);
    // The task-bar progress, created when the training starts: a request superseded while waiting shows none; the run's
    // record, built once the data is prepared, before the shared frame is given back.
    const run: {progress?: DG.TaskBarProgressIndicator; record?: TrainingRunRecord; fingerprint?: DatasetFingerprint;
      startedOn?: number} = {};
    try {
      this.updateControls();
      const queued = await this.queue.run(async (loop) => {
        run.startedOn = Date.now();
        const request = await prepareTraining(selection);
        try {
          run.fingerprint = datasetFingerprint(request.features, request.target);
          run.record = trainingRunOf({request, datasetName, fingerprint: run.fingerprint, status: 'completed',
            startedOn: new Date(run.startedOn).toISOString(), durationMs: 0});
          return await trainModel(request, loop);
        } finally {
          releaseFrame(request.features);
        }
      }, () => {
        run.progress = DG.TaskBarProgressIndicator.create(`Training ${engine.name} model`, {cancelable: true});
        return run.progress;
      });
      if (queued.outcome === 'superseded')
        return;
      const record = async (status: TrainingRunStatus, metrics?: MetricsRecord, message?: string) =>
        run.record === undefined ? '' : recordTrainingRun({...run.record, status, metrics, error: message,
          duration_ms: Date.now() - (run.startedOn ?? Date.now())});
      const isCurrent = () => checkNumber === this.checkNumber;
      if (queued.outcome === 'completed') {
        const {result} = queued;
        const runId = await record('completed', result.metrics);
        if (isCurrent() && run.fingerprint !== undefined) {
          this.lastTraining = {result, runId, engine, table, columns: [...features, target], datasetName,
            fingerprint: run.fingerprint, isSaved: false};
          this.isTrained = true;
          this.showResults(this.lastTraining);
        }
      } else if (queued.outcome === 'cancelled') {
        await record('cancelled');
        grok.shell.warning(CANCELLED);
        if (isCurrent())
          this.showHint(this.check?.isInteractive === false ? this.notInteractiveHint() : RESULTS_HINT);
      } else {
        await record('failed', undefined, errorMessage(queued.error));
        if (isCurrent()) {
          this.isTrained = true;
          this.showHint(RESULTS_HINT);
        }
        reportError(queued.error, isLive);
      }
    } catch (e) {
      reportError(e);
    } finally {
      run.progress?.close();
      this.activeTrainings--;
      if (!this.isTraining)
        ui.empty(this.resultsStatus);
      this.updateControls();
    }
  }

  private notInteractiveHint(): string {
    const rows = this.form.target.value?.length ?? this.form.table.rowCount;
    return `${this.engine.name} takes a while on ${rows} rows, so it does not retrain on every change. Click Train.`;
  }

  private showHint(text: string): void {
    ui.empty(this.resultsBody);
    this.resultsBody.append(ui.divText(text));
  }

  /** The loader next to the Results header; the previous results stay for a live training. */
  private showTraining(isLive: boolean): void {
    ui.empty(this.resultsStatus);
    this.resultsStatus.append(ui.loader());
    if (!isLive)
      ui.empty(this.resultsBody);
  }

  private showResults({result}: Training): void {
    ui.empty(this.resultsBody);
    this.resultsBody.append(metricsTable({metrics: result.metrics, rowCount: result.rowCount, folds: FOLDS,
      seed: result.seed, skippedRows: result.options.missingValues?.skippedRows ?? 0}));
  }

  private static valuesOf(inputs: ReadonlyMap<string, DG.InputBase>): Hyperparameters {
    const values: Hyperparameters = {};
    for (const [name, input] of inputs) {
      const value: unknown = input.value;
      if (typeof value === 'number' || typeof value === 'string' || typeof value === 'boolean')
        values[name] = value;
    }
    return values;
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

  /** The check's own problems: the data and the method, not the settings the user can fix in an input. */
  private static dataProblems({problems}: SelectionCheck): string[] {
    const {target, features, missingValues, method} = problems;
    return [...target, ...features, ...missingValues, ...method];
  }

  private static messageOf(messages: string[] | undefined): string | null {
    return messages?.length ? messages.join(' ') : null;
  }
}
