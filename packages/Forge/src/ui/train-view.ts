import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {ModelInfo} from '../catalog/model-edit';
import {defaultHyperparameters, enabledRules, EnabledRules, Engine, givesProbabilities, handlesMissingValues,
  Hyperparameters, isEnabled, isServerEngine, probabilityInput, visibleHyperparametersOf, withProbabilityInput}
  from '../engines/engine';
import {errorMessage} from '../forge-error';
import {ModelTask, TrainingRunStatus} from '../generated/db';
import {BestThresholdMetric, METRIC_LABELS, taskOf} from '../metrics/metrics';
import {_package} from '../package';
import {missingColumnsOf, MissingValuesSettings} from '../preparation/missing-values';
import {twoClasses} from '../preparation/pipeline';
import {DatasetFingerprint, datasetFingerprint, featuresHadGaps} from '../storage/dataset-fingerprint';
import {datasetRefOf} from '../storage/dataset-ref';
import {StorageChoice, TrainingRunRecord, trainingRunOf} from '../storage/model-fields';
import {saveTrainedModel} from '../storage/model-store';
import {recordTrainingRun} from '../storage/training-run-store';
import {browserSearch} from '../training/browser-search';
import {defaultFeatures} from '../training/default-features';
import {bestCandidates, bestSearchPlan, canContinue, failedSearchMessage, quickEngines, requestFingerprint,
  runPlanned, SearchChoice, SearchMode, searchableEngines, SearchState, searchTitle, tuneSearchPlan}
  from '../training/search';
import {bestThresholdOf, checkSelection, failedCheck, hasDataProblems, MetricsRecord, prepareTraining,
  problemMessages, recutResult, retrainsLive, selectedMissingValues, SelectionCheck, trainingEngines, TrainingRequest,
  TrainingResult, TrainingSelection, trainModel, trainPredictionOf} from '../training/train-model';
import {TrainingQueue} from '../training/training-queue';
import {openApplyDialog} from './apply-model-dialog';
import {ButtonGate} from './button-gate';
import {CollapsibleGroup} from './collapsible-group';
import {FindSettings} from './find-settings';
import {ForgeApp} from './forge-app';
import {MissingValuesInputs} from './missing-values-inputs';
import {ModelCards} from './model-cards';
import {PreparationInputs} from './preparation-inputs';
import {reportError} from './report-error';
import {ResultsPane} from './results-pane';
import {isSameInputs, RunInputs, RunsPanel, SessionRun} from './runs-panel';
import {saveModelDialog} from './save-model-dialog';
import {Restored, SearchPanel} from './search-panel';
import {MANUAL, TuneContext, TuneSettings, TuneStatus} from './tune-settings';
import {ValidationInputs} from './validation-inputs';

/** The built-in tool's delay after the last change before the preview is trained again. */
export const CHECK_DELAY_MS = 200;
const NAME_FEATURES_LENGTH = 30;
const RESULTS_HINT = 'Choose the target and the features; the model trains as you change them.';
const FIX_SETTINGS = 'Fix the settings.';
const FORM_WIDTH = 400;
const TRAIN_TOOLTIP = 'Train the model on the current selection.';
const TUNE_TOOLTIP = 'Find the best settings';
const FIND_TOOLTIP = 'Find the best model';
/** The first **Method** item: Find best model. */
const AUTO = 'Auto';
const METHOD_TOOLTIP = 'The machine learning method that builds the model.';
const AUTO_TOOLTIP = 'Auto: Find tries every checked method and trains the best one';
const AUTO_HINT = 'Click Find to try every checked method on this data and train the best one.';
const CONTINUE_TOOLTIP = 'Continue the search';
const CHECKING = 'Checking the selection...';
const TRAINING = 'Training is in progress.';
const TRAINED = 'The model is trained on this selection. Change an input to train again.';
const CANCELLED = 'Training was cancelled.';
const NO_TARGET = 'Choose a target.';
const TARGET_TOOLTIP = 'The column the model learns to predict.';
/** The id of the `training_run` row whose result the results show; absent while they show a hint. */
export const SHOWN_RUN = 'data-forge-run-id';

interface TrainForm {
  table: DG.DataFrame;
  tableInput: DG.InputBase<DG.DataFrame | null>;
  target: DG.InputBase<DG.Column | null>;
  features: DG.InputBase<DG.Column[]>;
  preparation: PreparationInputs;
  validation: ValidationInputs;
  method: DG.ChoiceInput<string | null>;
  hyperparameters: Map<string, DG.InputBase>;
  /** Settings (Manual / Auto-tune) and Search options of the method: rebuilt with the method. */
  tune: TuneSettings;
  /** The Find best model block: shown while Method is Auto. */
  find: FindSettings;
  /** The Method group's content: rebuilt with the method. */
  methodBody: HTMLDivElement;
  groups: CollapsibleGroup[];
}

export interface Training {
  result: TrainingResult;
  /** What the method got: the prepared columns, the selection's target and the rows kept. */
  request: TrainingRequest;
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
  /** **Train** right under the inputs, kept in view while the form scrolls: shown only for a method that does not
   * retrain live. */
  readonly trainButton: HTMLButtonElement;
  readonly saveButton: HTMLButtonElement;
  /** **Tune** under Search options in Auto-tune: searches the method's settings as they ask. */
  readonly tuneButton: HTMLButtonElement;
  /** **Find** under the Find best model block while Method is Auto: searches the checked methods. */
  readonly findButton: HTMLButtonElement;
  private readonly resultsStatus = ui.div([], 'forge-results-status');
  /** **Runs** under the results: the recorded trainings of this table in this view. */
  readonly runsPanel = new RunsPanel((run) => this.restoreRun(run));
  /** **Search** above **Runs**: the leaders of the last search of this table. */
  readonly searchPanel = new SearchPanel((leader) => this.restoreLeader(leader));
  /** The results: the **Performance** panel with the training loader in its header, the insights and the viewers; then
   * Search and Runs, all in one vertical splitter. */
  readonly resultsPane = new ResultsPane(this.resultsStatus, [{block: this.searchPanel.root, height: 320},
    {block: this.runsPanel.root, height: 260}]);
  /** The model cards under the form: the current model and the pinned ones (same data only). */
  readonly modelCards = new ModelCards((run) => this.restoreRun(run));
  lastTraining: Training | undefined;
  /** The recorded trainings of the current table, newest first; a Table change starts a new list. */
  runs: SessionRun[] = [];
  private readonly results: HTMLDivElement;
  private readonly trainRow: HTMLDivElement;
  private readonly trainGate: ButtonGate;
  private readonly tuneGate: ButtonGate;
  private readonly findGate: ButtonGate;
  /** The complete methods discovered when the view opened, in discovery order. */
  private readonly engines: Engine[];
  private readonly formPane: HTMLDivElement;
  private readonly changes = new rxjs.Subject<void>();
  private readonly queue = new TrainingQueue();
  /** Hyperparameter values per method name, kept for the session of the view (a Table change keeps them). */
  private readonly hyperparameterValues = new Map<string, Hyperparameters>();
  /** The values the last tune of each method set, the data it ran on and whether Find best model set them, kept as
   * {@link hyperparameterValues}: they read Auto-tuned on that data. */
  private readonly lastTuned = new Map<string, {values: Hyperparameters; context: TuneContext; isFound: boolean}>();
  /** The last search of each mode on this table (Tune's per method, {@link searchKey}) with the inputs it ran on
   * ({@link searchInputs}): **Tune** or **Find** continues it while they stay. In memory only: its plan holds the
   * columns. */
  private readonly searches = new Map<string, {state: SearchState; inputs: string}>();
  /** Methods whose check failed, logged once per view. */
  private readonly loggedFailures = new Set<string>();
  /** Whether each method retrains on every change, as the last check of this table that could ask said. */
  private readonly interactivity = new Map<string, boolean>();
  private engine: Engine;
  /** The method the user picked, or {@link AUTO}; null while the suggested one applies. */
  private userChoice: string | null = null;
  private form: TrainForm;
  // The current form's subscriptions: a Table change rebuilds the form and drops the old ones.
  private formSubs: rxjs.Subscription[] = [];
  // The hyperparameter inputs' subscriptions: a Method change rebuilds them.
  private methodSubs: rxjs.Subscription[] = [];
  /** The latest check's result; undefined while it is pending. */
  private check: SelectionCheck | undefined;
  private checkNumber = 0;
  /** An input other than the threshold changed since the timer last fired: the timer checks instead of re-cutting (as
   * it does when there is no result to re-cut and nothing trains). */
  private isCheckRequested = false;
  private activeTrainings = 0;
  /** A training of the current selection completed or failed; any change clears it. */
  private isTrained = false;
  private isFailed = false;
  /** Results being shown: `ResultsPane.show` calls not settled yet. */
  private showing = 0;
  private isSettingMethod = false;
  /** The current method's `meta.mlenabled` rules, read with the method. */
  private rules: EnabledRules = {};
  /** While a **Train** training or a search and its winner's training run: each form input's own enabled state, given
   * back at the end. */
  private locked: Map<DG.InputBase, boolean> | null = null;

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
    this.tuneButton = ui.button('Tune', () => void this.startSearch('tune', this.form.tune.choice()));
    this.tuneGate = new ButtonGate(this.tuneButton, () => this.tuneBlocker(), () => this.searchTooltip('tune'));
    this.findButton = ui.button('Find', () => void this.startSearch('best', this.form.find.choice()));
    this.findGate = new ButtonGate(this.findButton, () => this.findBlocker(), () => this.searchTooltip('best'));
    this.saveButton = ui.bigButton('Save', () => this.save());
    this.setRibbonPanels([[this.saveButton]]);
    this.formPane = ui.panel([], 'forge-train-form');
    this.results = ui.panel([this.resultsPane.root], 'forge-train-results');
    this.subs.push(DG.debounce(this.changes, CHECK_DELAY_MS).subscribe(() => {
      if (this.isCheckRequested || (this.lastTraining === undefined && !this.isTraining && !this.isTrained)) {
        this.isCheckRequested = false;
        void this.revalidate();
      } else
        this.recut();
    }));
    this.form = this.createForm(table);
    const left = ui.box(this.formPane);
    left.style.width = `${FORM_WIDTH}px`;
    this.root.appendChild(ui.splitH([left, ui.box(this.results)], {}, true));
  }

  getIcon(): HTMLElement {
    return ui.iconSvg('model');
  }

  detach(): void {
    this.queue.supersede();
    this.unsubscribeForm();
    this.unsubscribeMethod();
    this.resultsPane.clear();
    this.searchPanel.update(undefined);
    this.runsPanel.update([]);
    this.modelCards.clear();
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

  get tuneSettings(): TuneSettings {
    return this.form.tune;
  }

  /** The Find best model block of the form. */
  get findSettings(): FindSettings {
    return this.form.find;
  }

  /** **Method** is Auto: the Find best model block replaces the method's settings. */
  get isAutoMethod(): boolean {
    return this.userChoice === AUTO;
  }

  get missingValuesInputs(): MissingValuesInputs {
    return this.form.preparation.missingValues;
  }

  get preparationInputs(): PreparationInputs {
    return this.form.preparation;
  }

  get validationInputs(): ValidationInputs {
    return this.form.validation;
  }

  /** **Data**, **Preparation**, **Method** and **Validation**. */
  get groups(): readonly CollapsibleGroup[] {
    return this.form.groups;
  }

  /** A training is queued, running or being recorded. */
  get isTraining(): boolean {
    return this.activeTrainings > 0;
  }

  static async create(table: DG.DataFrame): Promise<TrainView> {
    const {engines, engine} = trainingEngines();
    return new TrainView(table, engines, engine);
  }

  /** Opens the view on the current table, or "Open a table first." without one: the menu's boundary. */
  static async open(): Promise<void> {
    const table = grok.shell.currentTable;
    if (table === null)
      grok.shell.warning('Open a table first.');
    else
      await TrainView.openFor(table);
  }

  /** Opens the view on [table]: the boundary of the menu and the table pane's **Train**. */
  static async openFor(table: DG.DataFrame): Promise<void> {
    try {
      grok.shell.addView(await TrainView.create(table));
    } catch (e) {
      reportError(e);
    }
  }

  /** **Train**: trains the current selection, whether the method retrains live or not; the form is locked meanwhile. */
  async train(): Promise<void> {
    if (this.problem() !== null || this.isTraining)
      return;
    this.lockForm();
    try {
      await this.startTraining(false);
    } finally {
      this.unlockForm();
    }
  }

  save(): void {
    const training = this.lastTraining;
    if (training === undefined || training.isSaved || this.form.preparation.thresholdProblem !== null)
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
    const {result, runId, datasetName, fingerprint, engine, table, columns} = training;
    let id: string;
    try {
      id = await saveTrainedModel({name, description, tags, engine, datasetName, result, fingerprint, columns, runId,
        storage: choice});
    } catch (e) {
      training.isSaved = false;
      this.updateControls();
      throw e;
    }
    grok.shell.info(ui.divV([ui.divText(`Model "${name}" saved.`), ui.divH([
      ui.link('Apply...', () => void openApplyDialog(table, {modelId: id, preferredTable: table})),
      ui.link('Show in the catalog', () => void ForgeApp.open()),
    ], 'forge-links')]));
  }

  /** The form for [table]; the Validation inputs keep the values of [previous] (the form of the table before). */
  private createForm(table: DG.DataFrame, previous?: ValidationInputs): TrainForm {
    const target = table.columns.byIndex(table.columns.length - 1);
    const onChanged = () => this.requestCheck();
    const preparation = new PreparationInputs(onChanged, () => this.onProbabilityToggled(), () => this.requestRecut(),
      (metric) => this.bestThreshold(metric), () => TrainView.messageOf(this.check?.problems.missingValues));
    const validation = new ValidationInputs(onChanged);
    if (previous !== undefined)
      validation.restore(previous.settings(), previous.seedValue);
    validation.update(target);
    const targetInput = ui.input.column('Target', {table, value: target, nullable: false,
      tooltipText: TrainView.targetTooltip(target), onValueChanged: (t, input) => {
        input.setTooltip(TrainView.targetTooltip(t));
        preparation.update(featuresInput.value, t);
        validation.update(t);
        onChanged();
      }});
    const checked = defaultFeatures(table, target).map((c) => c.name);
    const featuresInput = ui.input.columns('Features', {table, checked, nullable: false,
      tooltipText: 'Columns the model uses to make predictions.', onValueChanged: (columns) => {
        preparation.update(columns, targetInput.value);
        onChanged();
      }});
    preparation.update(featuresInput.value, target);
    targetInput.addValidator(() => TrainView.messageOf(this.check?.problems.target));
    featuresInput.addValidator(() => TrainView.messageOf(this.check?.problems.features));
    const tableInput = ui.input.table('Table', {items: grok.shell.tables, value: table,
      tooltipText: 'Data to train the model on.', onValueChanged: (t) => {
        if (t !== null)
          this.form = this.createForm(t, this.form.validation);
      }});
    const methodInput = ui.input.choice<string | null>('Method', {items: [AUTO, this.engine.name],
      value: this.isAutoMethod ? AUTO : this.engine.name, nullable: false, tooltipText: METHOD_TOOLTIP,
      onValueChanged: (name) => this.chooseMethod(name)});
    methodInput.addValidator(() => TrainView.messageOf(this.check?.problems.method));
    const find = new FindSettings(this.findButton, methodInput.captionLabel, {onChanged: () => {
      if (this.isAutoMethod)
        preparation.setProbabilityReason(find.probabilityReason);
      this.updateControls();
    }});

    const dataInputs = [tableInput, targetInput, featuresInput];
    const methodBody = ui.div([]);
    const groups = [new CollapsibleGroup('Data', dataInputs.map((input) => input.root)), preparation.group,
      new CollapsibleGroup('Method', [methodBody]), validation.group];
    this.unsubscribeForm();
    this.interactivity.clear();
    ui.setDisplay(this.trainRow, false);
    this.formSubs = [...groups[0].expandOnError(dataInputs), ...groups[1].expandOnError(preparation.inputs),
      ...groups[2].expandOnError([methodInput, ...find.inputs]), ...groups[3].expandOnError(validation.inputs),
      ...find.subs,
      // The platform's column inputs keep a removed column until their dialog is opened again.
      table.onColumnsRemoved.subscribe(() => TrainView.dropRemovedColumns(table, targetInput, featuresInput)),
      // The task the method trains decides some of its settings' rules.
      targetInput.onChanged.subscribe(() => this.applyRules(this.form)),
      rxjs.merge(targetInput.onChanged, featuresInput.onChanged).subscribe(() => this.form.tune.updateStatus())];
    this.modelCards.clear();
    ui.empty(this.formPane);
    // One form, so the labels of every group line up.
    const formRoot = ui.form([]);
    formRoot.append(...groups.map((g) => g.root), this.trainRow, this.modelCards.root);
    this.formPane.append(formRoot);
    this.lastTraining = undefined;
    this.runs = [];
    this.searches.clear();
    this.runsPanel.update(this.runs);
    this.searchPanel.update(undefined);
    this.resultsPane.clear();
    this.resultsPane.layout();
    this.showHint(RESULTS_HINT);
    // fillMethod adds the tune settings.
    const form = {table, tableInput, target: targetInput, features: featuresInput, preparation, validation,
      method: methodInput, hyperparameters: new Map(), find, methodBody, groups} as TrainForm;
    this.fillMethod(form);
    this.requestCheck();
    return form;
  }

  /** The Method group of [form] for the current method: Method, its Settings (Manual / Auto-tune) and its
   * hyperparameters with the values kept for it. */
  private fillMethod(form: TrainForm): void {
    const engine = this.engine;
    const rules = enabledRules(engine);
    const key = `${engine.name} mlenabled`;
    if (rules === undefined && !this.loggedFailures.has(key)) {
      this.loggedFailures.add(key);
      _package.logger.warning(`The method ${engine.name} declares a malformed meta.mlenabled; ` +
        'every setting stays enabled.');
    }
    this.rules = rules ?? {};
    // Under Auto the search never keeps the rows with missing values, whatever the suggested method could.
    form.preparation.setKeepsMissingValues(!this.isAutoMethod && handlesMissingValues(engine));
    const {predictProbability} = form.preparation;
    form.find.setPredictProbability(predictProbability.value);
    form.preparation.setProbabilityReason(this.probabilityReason(form));
    form.method.setTooltip(this.isAutoMethod ? AUTO_TOOLTIP : METHOD_TOOLTIP);
    if (this.locked !== null) {
      this.locked.set(predictProbability, predictProbability.enabled);
      predictProbability.enabled = false;
    }
    const values = {...defaultHyperparameters(engine), ...this.hyperparameterValues.get(engine.name)};
    this.unsubscribeMethod();
    form.hyperparameters.clear();
    for (const p of visibleHyperparametersOf(engine)) {
      const input = ui.input.forProperty(p);
      const value = values[p.name];
      if (value !== undefined)
        input.value = value;
      form.hyperparameters.set(p.name, input);
      this.methodSubs.push(input.onChanged.subscribe(() => {
        this.hyperparameterValues.set(engine.name, TrainView.valuesOf(form.hyperparameters));
        this.applyRules(form);
        this.requestCheck();
      }));
    }
    form.tune = new TuneSettings(engine, form.hyperparameters, this.tuneButton, {
      values: () => TrainView.valuesOf(form.hyperparameters),
      context: () => form.target.value === null ? undefined : {table: form.table, target: form.target.value.name,
        features: form.features.value.map((c) => c.name)},
      lastTuned: () => this.lastTuned.get(engine.name),
      onChanged: (affectsRules) => {
        if (affectsRules)
          this.applyRules(form);
        this.showTrainRow();
        this.updateControls();
      }});
    const inputs = [...form.hyperparameters.values()];
    this.methodSubs.push(...form.tune.subs, ...form.groups[2].expandOnError([...inputs, ...form.tune.inputs]));
    ui.empty(form.methodBody);
    form.methodBody.append(form.method.root, form.tune.root, form.find.root);
    ui.setDisplay(form.tune.root, !this.isAutoMethod);
    ui.setDisplay(form.find.root, this.isAutoMethod);
    if (this.locked !== null) {
      form.tune.setLocked(true);
      form.find.setLocked(true);
    }
    this.applyRules(form);
  }

  /** Why Predict probability cannot be used: under Auto no checked method gives class probabilities, else the method
   * gives none; null when it can. */
  private probabilityReason(form: TrainForm): string | null {
    return this.isAutoMethod ? form.find.probabilityReason : givesProbabilities(this.engine) ? null :
      `${this.engine.name} gives no class probabilities`;
  }

  /** Disables the settings of [form] the method ignores with its values and task (`meta.mlenabled`) and enables the
   * others; a disabled setting keeps its value. In Auto-tune a rule on a tuned setting may hold. The setting's row in
   * Search options follows it. While the form is locked, the states are kept for the unlock. */
  private applyRules(form: TrainForm): void {
    const {tune} = form;
    const target = form.target.value;
    const task: ModelTask | undefined = target === null ? undefined : taskOf(target);
    tune.setTask(task);
    form.find.setTask(task);
    const free = tune.isAuto ? tune.tunedNames : undefined;
    const values = TrainView.valuesOf(form.hyperparameters);
    const inactive = new Set<string>();
    for (const [name, input] of form.hyperparameters) {
      const isActive = isEnabled(this.rules, name, values, task, free);
      if (!isActive)
        inactive.add(name);
      this.locked?.set(input, isActive);
      const enabled = isActive && this.locked === null;
      if (input.enabled !== enabled)
        input.enabled = enabled;
    }
    tune.setInactive(inactive);
  }

  /** The user picked [name] in **Method**: a method, or Auto, which shows the Find best model block and keeps the
   * shown model. */
  private chooseMethod(name: string | null): void {
    if (this.isSettingMethod)
      return;
    const wasAuto = this.isAutoMethod;
    if (name === AUTO) {
      this.userChoice = AUTO;
      if (!wasAuto)
        this.fillMethod(this.form);
      this.requestCheck(true);
      return;
    }
    const engine = this.engines.find((e) => e.name === name);
    if (engine === undefined)
      return;
    this.userChoice = engine.name;
    if (engine !== this.engine || wasAuto) {
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

  private selection(): TrainingSelection | null {
    const target = this.form.target.value;
    if (target === null)
      return null;
    const {validation} = this.form;
    const steps = this.form.preparation.steps();
    return {engine: this.engine, features: this.form.features.value, target,
      hyperparameters: withProbabilityInput(this.engine, TrainView.valuesOf(this.form.hyperparameters),
        steps.predictProbability), seed: validation.seedValue,
      validation: validation.settings(), missingValues: this.form.preparation.missingValues.settings(), steps};
  }

  /** Predict probability changed: under Auto the Candidates list follows it and nothing trains; else a method whose
   * probabilities are a setting trains again, any other re-cuts the result. */
  private onProbabilityToggled(): void {
    if (this.isAutoMethod) {
      this.form.find.setPredictProbability(this.form.preparation.steps().predictProbability);
      this.updateControls();
    } else if (probabilityInput(this.engine) !== undefined)
      this.requestCheck();
    else
      this.requestRecut();
  }

  /** A change that keeps the result: the controls follow it, and the timer re-cuts the result. */
  private requestRecut(): void {
    this.updateControls();
    this.changes.next();
  }

  /** The threshold that maximizes [metric] on the validation rows of the last training ({@link bestThresholdOf});
   * undefined without a training, validation probabilities or two classes in the validation rows. */
  private bestThreshold(metric: BestThresholdMetric): number | undefined {
    const result = this.lastTraining?.result;
    return result === undefined ? undefined : bestThresholdOf(result, metric);
  }

  /** Every change: the training of the old selection stops, its result goes unless [keepsResult] (choosing Auto keeps
   * the selection), and the selection is checked again [CHECK_DELAY_MS] after the last change. */
  private requestCheck(keepsResult = false): void {
    this.checkNumber++;
    this.isCheckRequested = true;
    this.check = undefined;
    this.queue.supersede();
    if (!keepsResult) {
      this.lastTraining = undefined;
      this.isTrained = false;
      this.isFailed = false;
    }
    this.updateControls();
    this.changes.next();
  }

  private async checkCurrent(): Promise<SelectionCheck> {
    const selection = this.selection();
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
      if (engine !== this.engine && check.prepared !== undefined) {
        this.engine = engine;
        this.fillMethod(this.form);
        // The data and the list stand; only the new method, which is listed, is asked whether it retrains live.
        const {columns, target} = check.prepared;
        const isInteractive = await retrainsLive(engine, columns, target);
        if (checkNumber !== this.checkNumber)
          return;
        check = {...check, problems: {...check.problems, method: []}, isInteractive};
      }
      // Under Auto the suggested method is not trained: its "cannot learn" is no problem while another method can.
      if (this.isAutoMethod && check.engines.length > 0)
        check = {...check, problems: {...check.problems, method: []}};
      if (this.isAutoMethod && check.prepared !== undefined) {
        const {columns, target} = check.prepared;
        const browser = searchableEngines(check.engines, this.engines);
        const live = await quickEngines(browser, columns, target, false);
        if (checkNumber !== this.checkNumber)
          return;
        for (const e of browser)
          this.interactivity.set(e.name, live.includes(e));
        this.form.find.setMethods([...browser.map((e) => FindSettings.entryOf(e, live.includes(e))),
          ...check.engines.filter(isServerEngine).map((e) => FindSettings.entryOf(e, false))]);
        this.form.preparation.setProbabilityReason(this.probabilityReason(this.form));
      }
      // The prepared columns (one-hot copies of the whole table) are not kept beyond this check.
      this.check = {...check, prepared: undefined};
      this.resultsPane.setDataInsights(check.insights);
      this.showMethods(check);
      const settings = this.form.preparation.missingValues.visibleInputs;
      for (const input of [this.form.target, this.form.features, this.form.method, this.form.preparation.threshold,
        ...settings, ...this.form.validation.usedInputs, ...this.activeHyperparameters(), ...this.form.find.inputs])
        input.validate();
      if (check.engines.length > 0)
        this.interactivity.set(this.engine.name, check.isInteractive);
      this.showTrainRow();
      this.updateControls();
      if (this.problem() !== null)
        this.showHint(TrainView.dataProblems(check).length === 0 ? FIX_SETTINGS : RESULTS_HINT);
      else if (this.isAutoMethod) {
        if (this.lastTraining === undefined)
          this.showHint(AUTO_HINT);
      } else if (check.isInteractive)
        await this.startTraining(true);
      else
        this.showHint(this.notInteractiveHint());
    } catch (e) {
      if (checkNumber !== this.checkNumber)
        return;
      this.check = failedCheck({method: [errorMessage(e)]});
      this.resultsPane.setDataInsights([]);
      this.showTrainRow();
      this.updateControls();
      reportError(e, true);
    }
  }

  /** **Train** is shown in Manual for a method known not to retrain on every change, never for an interactive or an
   * unknown one, nor under Auto: the data checks that fail leave the method's last interactivity on this table. */
  private showTrainRow(): void {
    ui.setDisplay(this.trainRow, this.interactivity.get(this.engine.name) === false && !this.form.tune.isAuto &&
      !this.isAutoMethod);
  }

  /** The method for [check]: the current one while a target or features rule fails or nothing can learn; the
   * suggested one under Auto; the user's choice while it is listed; the suggested one otherwise, with a balloon when it
   * replaces the user's choice. */
  private methodOf(check: SelectionCheck): Engine {
    const {problems, engines, best} = check;
    if (hasDataProblems(problems) || engines.length === 0)
      return this.engine;
    if (this.isAutoMethod)
      return best ?? engines[0];
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
      this.form.method.items = [AUTO, ...engines.map((e) => e.name)];
      this.form.method.value = this.isAutoMethod ? AUTO : engines.includes(this.engine) ? this.engine.name : null;
    } finally {
      this.isSettingMethod = false;
    }
  }

  /** Why **Train** is unavailable: the check is pending, a problem, a training running or the selection already
   * trained; null when it is not. */
  private blocker(problem: string | null = this.problem()): string | null {
    return problem ?? (this.isTraining ? TRAINING : this.isTrained ? TRAINED : null);
  }

  /** Why **Tune** is unavailable: as **Train**, but a trained selection can be tuned; a method without a search
   * space cannot, nor Search options with a problem. */
  private tuneBlocker(problem: string | null = this.problem()): string | null {
    return problem ?? (this.isTraining ? TRAINING : this.form.tune.noSpaceReason ?? this.form.tune.problem());
  }

  /** Why **Find** is unavailable: as **Tune**, or the Find best model block's problem. */
  private findBlocker(problem: string | null = this.problem()): string | null {
    return problem ?? (this.isTraining ? TRAINING : this.form.find.problem());
  }

  /** The pending check or the selection's first problem, an invalid validation input or hyperparameter included (one
   * the method ignores is not checked; none under Auto); null when the selection can be trained. */
  private problem(): string | null {
    if (this.check === undefined)
      return CHECKING;
    const captioned = (inputs: DG.InputBase[]) => inputs.filter((input) => input.validity !== null)
      .map((input) => `${input.caption}: ${input.validity}`);
    const settings = this.form.preparation.missingValues.visibleInputs.map((input) => input.validity)
      .filter((v): v is string => v !== null);
    const threshold = this.form.preparation.thresholdProblem;
    return [...TrainView.dataProblems(this.check), ...settings, ...captioned(this.form.validation.usedInputs),
      ...(this.isAutoMethod ? [] : captioned(this.activeHyperparameters())),
      ...(threshold === null ? [] : [`Classification threshold: ${threshold}`])][0] ?? null;
  }

  /** The hyperparameter inputs of the settings the method uses with the form's values. */
  private activeHyperparameters(): DG.InputBase[] {
    return [...this.form.hyperparameters].filter(([name]) => this.form.tune.isActive(name)).map(([, input]) => input);
  }

  /** **Train**, **Tune** (in Auto-tune only: a hidden one keeps no disabled state), **Find** (under Auto only),
   * **Save** and the Results hint follow the check, the training and the last result. */
  private updateControls(): void {
    const problem = this.problem();
    this.trainGate.update(this.blocker(problem));
    // The first form is still being built while its check is requested: there is no form yet, and no Auto.
    const isTune = !this.isAutoMethod && this.form?.tune.isAuto === true;
    if (this.isAutoMethod)
      this.findGate.update(this.findBlocker(problem));
    else
      this.findGate.release();
    if (isTune)
      this.tuneGate.update(this.tuneBlocker(problem));
    else
      this.tuneGate.release();
    const target = this.form?.target.value;
    this.modelCards.update(target == null ? undefined :
      ModelCards.keyOf(target.name, this.form.validation.settings(), this.form.validation.seedValue));
    const training = this.lastTraining;
    ui.setDisabled(this.saveButton, this.isTraining || training === undefined || training.isSaved ||
      this.form.preparation.thresholdProblem !== null);
    this.form?.preparation.updateSet();
    this.publishState(problem);
  }

  /** The results' readiness for automation (ARCHITECTURE.md, "Automation surface"). */
  private publishState(problem: string | null): void {
    TrainView.setAttribute(this.results, 'aria-busy', String(this.check === undefined || this.isTraining ||
      this.showing > 0));
    TrainView.setAttribute(this.results, 'aria-invalid', String(this.check !== undefined && (problem !== null ||
      this.isFailed)));
  }

  private static setAttribute(element: Element, name: string, value: string): void {
    if (element.getAttribute(name) !== value)
      element.setAttribute(name, value);
  }

  /** Trains the current selection through the queue (a newer training supersedes it) with a task-bar progress; records
   * a completed, failed or cancelled run (a superseded one is not recorded) and shows the result. [override] is the
   * winner of a search: its selection, its request (the search rows, the hold out and the search record) and the
   * fingerprint of all the selection's rows. */
  private async startTraining(isLive: boolean,
    override?: {selection: TrainingSelection; request: TrainingRequest; fingerprint: DatasetFingerprint}):
    Promise<void> {
    const selection = override?.selection ?? this.selection();
    if (selection === null)
      return;
    const settingsStatus = this.form.tune.status;
    const checkNumber = this.checkNumber;
    // A Table change starts a new list: a training of the old table adds no row to it.
    const runs = this.runs;
    const {engine, features, target} = selection;
    const table = this.form.table;
    const datasetName = table.name;
    this.activeTrainings++;
    this.showTraining(isLive);
    // The task-bar progress, created when the training starts: a request superseded while waiting shows none; the run's
    // record, built once the data is prepared, before the shared frame is given back.
    const run: {progress?: DG.TaskBarProgressIndicator; request?: TrainingRequest; record?: TrainingRunRecord;
      fingerprint?: DatasetFingerprint; startedOn?: number} = {};
    try {
      this.updateControls();
      const queued = await this.queue.run(async (loop) => {
        run.startedOn = Date.now();
        const request = override?.request ?? await prepareTraining(selection);
        run.request = request;
        run.fingerprint = override?.fingerprint ?? datasetFingerprint(request.features, request.target);
        run.record = trainingRunOf({request, datasetName, fingerprint: run.fingerprint, status: 'completed',
          source: 'ui', startedOn: new Date(run.startedOn).toISOString(), durationMs: 0});
        return await trainModel(request, loop);
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
      const addRun = (status: TrainingRunStatus, outcome: {metrics?: MetricsRecord; error?: string},
        steps = selection.steps) => {
        if (runs === this.runs)
          this.addRun({...selection, steps}, run.record, status, outcome, settingsStatus);
      };
      if (queued.outcome === 'completed') {
        const {result} = queued;
        const {binaryClassificationThreshold, thresholdMetric} = result.options;
        if (run.record !== undefined)
          run.record.options = result.options;
        const runId = await record('completed', result.metrics);
        addRun('completed', {metrics: result.metrics}, {...selection.steps,
          threshold: binaryClassificationThreshold ?? selection.steps.threshold, thresholdMetric});
        if (isCurrent() && run.fingerprint !== undefined && run.request !== undefined) {
          this.lastTraining = {result, request: run.request, runId, engine, table, columns: [...features, target],
            datasetName, fingerprint: run.fingerprint, isSaved: false};
          this.isTrained = true;
          this.isFailed = false;
          // Predict probability or the threshold may have changed while the model trained.
          this.recut();
        }
      } else if (queued.outcome === 'cancelled') {
        await record('cancelled');
        addRun('cancelled', {});
        grok.shell.warning(CANCELLED);
        if (isCurrent()) {
          this.isFailed = true;
          this.idleHint();
        }
      } else {
        const error = errorMessage(queued.error);
        await record('failed', undefined, error);
        addRun('failed', {error});
        if (isCurrent()) {
          this.isTrained = true;
          this.isFailed = true;
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

  /** Puts the training of [selection] whose outcome was recorded as [record] (undefined when nothing was) first in
   * **Runs**, numbered after the runs before it; a completed one becomes the current model card. */
  private addRun(selection: TrainingSelection, record: TrainingRunRecord | undefined, status: TrainingRunStatus,
    outcome: {metrics?: MetricsRecord; error?: string}, settingsStatus: TuneStatus): void {
    this.runs.unshift({number: this.runs.length + 1, status, ...TrainView.runInputsOf(selection),
      featureGaps: missingColumnsOf(selection.features).length > 0, rowCount: record?.row_count ?? 0,
      search: record?.search, settingsStatus, ...outcome});
    this.runsPanel.update(this.runs);
    this.resultsPane.layout();
    if (status === 'completed')
      this.modelCards.setCurrent(this.runs[0]);
  }

  /** Sets the target, the features, the preparation, the validation, the method and its hyperparameters of [run]
   * again; the view then checks and trains as after any change. A run whose columns the table no longer has sets
   * nothing with a warning; a run whose inputs the form has and shows the completed result of sets nothing, or only
   * its Predict probability and threshold, which re-cut that result. The settings switch to Manual. Nothing is set
   * while the form is locked. */
  restoreRun(run: SessionRun): void {
    if (this.locked !== null)
      return;
    const table = this.form.table;
    const target = table.col(run.targetName);
    const missing = [run.targetName, ...run.featureNames].filter((name) => table.col(name) === null);
    if (target === null || missing.length > 0) {
      grok.shell.warning(`The table no longer has the columns of this run: ${missing.join(', ')}.`);
      return;
    }
    this.form.tune.setMode(MANUAL);
    // Any change but Predict probability (of a method without a probability setting) and the threshold drops
    // lastTraining: while it is there, it is the result of the form's inputs.
    const current = this.selection();
    if (this.lastTraining !== undefined && current !== null && isSameInputs(run, TrainView.runInputsOf(current))) {
      const {preparation} = this.form;
      if (preparation.predictProbability.value !== run.steps.predictProbability)
        preparation.predictProbability.value = run.steps.predictProbability;
      preparation.restoreThreshold(run.steps);
      this.requestRecut();
      return;
    }
    this.setMethod(run.engine, run.hyperparameters);
    this.form.target.value = target;
    this.form.features.value = run.featureNames.map((name) => table.getCol(name));
    this.form.preparation.restore(run.steps, run.missingValues);
    this.form.validation.restore(run.validation, run.seed);
    this.requestCheck();
  }

  /** Searches on the selection as [choice] asks, through the queue with a task-bar progress: [mode] `tune` the current
   * method's settings (**Tune**: as Search options set it), `best` the methods [choice] names (**Find**: as the Find
   * best model block sets it; the rows with missing values skipped unless imputed). Shows the leaders in **Search**,
   * sets the form to the chosen method and settings, which read Auto-tuned in Manual, and trains them as a regular
   * training on the search rows, with the hold out measured as Test. The kept search of [mode] on the same inputs
   * ({@link keptSearch}) is continued when its plan allows ({@link canContinue}): the new budget on top of its trials;
   * otherwise the search starts afresh; either way it becomes the kept one. Meanwhile the form is locked: its inputs
   * are disabled, and a click on a **Runs** or **Search** row sets nothing. */
  async startSearch(mode: SearchMode, choice: SearchChoice): Promise<void> {
    const selection = this.searchSelection(mode);
    if (selection === null || (mode === 'tune' ? this.tuneBlocker() : this.findBlocker()) !== null)
      return;
    const {engine} = this;
    const {space} = this.form.tune;
    const {seed} = selection;
    const key = this.searchKey(mode);
    const inputs = this.searchInputs(mode, selection, choice);
    const previous = this.keptSearch(mode, selection, choice);
    const title = searchTitle(mode, engine.name);
    const checkNumber = this.checkNumber;
    this.activeTrainings++;
    this.showTraining(false);
    this.searchPanel.update(undefined);
    this.resultsPane.layout();
    this.lockForm();
    let indicator: DG.TaskBarProgressIndicator | undefined;
    try {
      this.updateControls();
      let queued;
      try {
        queued = await this.queue.run(async (loop) => {
          const request = await prepareTraining(selection);
          const fingerprint = datasetFingerprint(request.features, request.target);
          const planned = mode === 'tune' ?
            tuneSearchPlan(request, engine, selection.hyperparameters, choice, seed, space) :
            bestSearchPlan(request, bestCandidates(this.candidateEngines(choice.methods ?? []), this.engines,
              selection.steps.predictProbability), choice, seed);
          const requestKey = requestFingerprint(planned.plan.request);
          const continued = previous !== undefined && canContinue(previous, planned.plan, requestKey) ? previous :
            undefined;
          if (continued !== undefined)
            loop.update(0, `Continuing: ${title}`);
          return {fingerprint, ...await runPlanned(planned, loop, browserSearch, continued, requestKey)};
        }, () => {
          indicator = DG.TaskBarProgressIndicator.create(title, {cancelable: true});
          return indicator;
        });
      } finally {
        indicator?.close();
      }
      if (queued.outcome === 'superseded')
        return;
      const isCurrent = checkNumber === this.checkNumber;
      if (queued.outcome !== 'completed') {
        if (isCurrent)
          this.idleHint();
        reportError(queued.error);
        return;
      }
      const {plan, holdOut, outcome, leaders, chosen, record, fingerprint, state} = queued.result;
      if (!isCurrent)
        return;
      this.searches.set(key, {state, inputs});
      this.searchPanel.update({plan, leaders, outcome});
      this.resultsPane.layout();
      if (chosen === undefined) {
        grok.shell.warning(failedSearchMessage(outcome) ?? 'No trial completed.');
        this.idleHint();
        return;
      }
      const winner = plan.candidates.find((c) => c.engine.name === chosen.method)!.engine;
      this.setWinner(winner, chosen.hyperparameters);
      this.lastTuned.set(winner.name, {values: TrainView.valuesOf(this.form.hyperparameters), context: {
        table: this.form.table, target: selection.target.name, features: selection.features.map((c) => c.name)},
      isFound: mode === 'best'});
      this.form.tune.updateStatus();
      if (outcome.stopped === 'cancelled')
        grok.shell.warning('Search cancelled; the best so far is trained.');
      else if (mode === 'tune')
        grok.shell.info(`Optimal settings for ${chosen.method} found.`);
      else
        grok.shell.info(`Best model: ${chosen.method} (${METRIC_LABELS[plan.metric]} ${chosen.score.toFixed(3)}).`);
      const hyperparameters = {...chosen.hyperparameters};
      await this.startTraining(false, {selection: {...selection, engine: winner, hyperparameters},
        request: {...plan.request, engine: winner, hyperparameters, holdOut, search: record}, fingerprint});
    } catch (e) {
      reportError(e);
    } finally {
      this.unlockForm();
      this.activeTrainings--;
      if (!this.isTraining)
        ui.empty(this.resultsStatus);
      this.updateControls();
    }
  }

  /** The selection a search of [mode] runs on: Find best model's with its missing values ({@link searchMissingValues});
   * null without a target. */
  private searchSelection(mode: SearchMode): TrainingSelection | null {
    const base = this.selection();
    return base === null || mode === 'tune' ? base : {...base, missingValues: TrainView.searchMissingValues(base)};
  }

  /** The key of the kept search of [mode]: Tune's per method. */
  private searchKey(mode: SearchMode): string {
    return mode === 'tune' ? `tune ${this.engine.name}` : 'best';
  }

  /** What a search of [mode] on [selection] as [choice] asks continues only while it stays the same: the target, the
   * features, the preparation (the threshold too, which the trials score with: {@link requestFingerprint}), the missing
   * values, the validation and its seed, the metric, the strategy and the methods (Find best model), the space and the
   * values of the settings not tuned (Tune). The limit is not part of it: a continuation adds it. */
  private searchInputs(mode: SearchMode, selection: TrainingSelection, choice: SearchChoice): string {
    const {target, features, hyperparameters, steps, validation, seed} = selection;
    const tuned = this.form.tune.tunedNames;
    const fixed = mode === 'tune' ? Object.entries(hyperparameters).filter(([name]) => !tuned.has(name)) : [];
    return JSON.stringify([target.name, features.map((c) => c.name), steps, selectedMissingValues(selection),
      validation, seed, choice.metric ?? null, choice.strategy ?? null, choice.methods ?? [], choice.space ?? null,
      fixed]);
  }

  /** The kept search of [mode] when [selection] and [choice] ask what it ran on ({@link searchInputs}). */
  private keptSearch(mode: SearchMode, selection: TrainingSelection, choice: SearchChoice): SearchState | undefined {
    const kept = this.searches.get(this.searchKey(mode));
    return kept !== undefined && kept.inputs === this.searchInputs(mode, selection, choice) ? kept.state : undefined;
  }

  /** The tooltip of **Tune** ([mode] `tune`) or **Find**: the trials so far of the search it would continue, else what
   * it does. */
  private searchTooltip(mode: SearchMode): string {
    const selection = this.searches.size === 0 ? null : this.searchSelection(mode);
    const kept = selection === null ? undefined :
      this.keptSearch(mode, selection, mode === 'tune' ? this.form.tune.choice() : this.form.find.choice());
    return kept !== undefined ? `${CONTINUE_TOOLTIP}: ${kept.outcome.trials.length} trials so far` :
      mode === 'tune' ? TUNE_TOOLTIP : FIND_TOOLTIP;
  }

  /** The browser methods of the last check named in [names], cheapest first. */
  private candidateEngines(names: string[]): Engine[] {
    return searchableEngines(this.check?.engines ?? [], this.engines).filter((e) => names.includes(e.name));
  }

  /** Find best model's missing values: imputed when the form imputes, else skipped; never kept, as every method
   * gets the same folds. */
  private static searchMissingValues(selection: TrainingSelection): MissingValuesSettings {
    const chosen = selectedMissingValues(selection);
    return chosen.mode === 'impute' ? chosen : {mode: 'skip'};
  }

  /** Every input of the form: the data, the preparation, the validation, the method, its settings and the Find best
   * model block. */
  private formInputs(): DG.InputBase[] {
    const {tableInput, target, features, preparation, validation, method, hyperparameters, find} = this.form;
    return [tableInput, target, features, ...preparation.inputs, ...validation.inputs, method,
      ...hyperparameters.values(), ...find.inputs];
  }

  /** Disables every input of the form, Settings and the controls of Search options, keeping each input's own enabled
   * state for {@link unlockForm}. Locks never nest: both callers start only while nothing trains. */
  private lockForm(): void {
    const locked = new Map<DG.InputBase, boolean>();
    for (const input of this.formInputs()) {
      locked.set(input, input.enabled);
      input.enabled = false;
    }
    this.locked = locked;
    this.form.tune.setLocked(true);
    this.form.find.setLocked(true);
  }

  /** Gives each input of the form its enabled state back, the method's settings by their rules. */
  private unlockForm(): void {
    const locked = this.locked;
    this.locked = null;
    for (const [input, enabled] of locked ?? [])
      input.enabled = enabled;
    this.form.tune.setLocked(false);
    this.form.find.setLocked(false);
    this.applyRules(this.form);
  }

  /** Sets **Method** to [engine] with [hyperparameters], without a check: the search's winner is trained next. */
  private setWinner(engine: Engine, hyperparameters: Hyperparameters): void {
    this.setMethod(engine, hyperparameters);
    this.isSettingMethod = true;
    try {
      if (!this.form.method.items.includes(engine.name))
        this.form.method.items = [...this.form.method.items, engine.name];
      this.form.method.value = engine.name;
    } finally {
      this.isSettingMethod = false;
    }
    this.showTrainRow();
  }

  /** A click on a **Search** row: sets its method and settings in Manual; the view then checks and trains as after
   * any change. Nothing is set while the form is locked. */
  restoreLeader(leader: Restored): void {
    const engine = this.engines.find((e) => e.name === leader.method);
    if (engine === undefined || this.locked !== null)
      return;
    this.setMethod(engine, leader.hyperparameters);
    this.requestCheck();
  }

  /** Sets **Method** to [engine] with [hyperparameters] in Manual. */
  private setMethod(engine: Engine, hyperparameters: Hyperparameters): void {
    this.hyperparameterValues.set(engine.name, {...hyperparameters});
    this.userChoice = engine.name;
    this.engine = engine;
    this.fillMethod(this.form);
    this.form.tune.setMode(MANUAL);
  }

  /** The Results hint while nothing is shown: under Auto the last training when there is one, else Click Find; Click
   * Train for a method that does not retrain live. */
  private idleHint(): void {
    if (this.isAutoMethod && this.lastTraining !== undefined)
      this.showResults(this.lastTraining);
    else {
      this.showHint(this.isAutoMethod ? AUTO_HINT : this.check?.isInteractive === false ? this.notInteractiveHint() :
        RESULTS_HINT);
    }
  }

  private notInteractiveHint(): string {
    const rows = this.form.target.value?.length ?? this.form.table.rowCount;
    return `${this.engine.name} takes a while on ${rows} rows, so it does not retrain on every change. Click Train.`;
  }

  private showHint(text: string): void {
    this.resultsPane.showHint(text);
    this.results.removeAttribute(SHOWN_RUN);
  }

  /** The loader next to the Performance header; the previous results stay for a live training. */
  private showTraining(isLive: boolean): void {
    ui.empty(this.resultsStatus);
    this.resultsStatus.append(ui.loader());
    if (!isLive)
      this.showHint('');
  }

  /** Shows the last training with the current **Predict probability** and **Classification threshold**, without
   * retraining; Save writes it so. A threshold not edited by hand becomes the chosen metric's best on this result. */
  private recut(): void {
    const training = this.lastTraining;
    if (training === undefined)
      return;
    const {preparation} = this.form;
    preparation.threshold.validate();
    if (preparation.thresholdProblem !== null) {
      this.showHint(FIX_SETTINGS);
      return;
    }
    training.result = recutResult(training.result, preparation.steps());
    const {binaryClassificationThreshold, thresholdMetric} = training.result.options;
    if (binaryClassificationThreshold !== undefined && thresholdMetric !== undefined)
      preparation.setThreshold(binaryClassificationThreshold, thresholdMetric);
    this.updateControls();
    this.showResults(training);
  }

  /** Shows [training] in **Results**: its metrics (in the grid it already holds while the metric rows stay), the
   * insights and the viewers over the rows it was trained on. */
  private showResults({result, request, engine, table, runId, fingerprint}: Training): void {
    const summary = {metrics: result.metrics, rowCount: result.rowCount, splitting: result.splitting, seed: result.seed,
      skippedRows: result.options.missingValues?.skippedRows ?? 0,
      keepsMissingValues: result.options.missingValues?.mode === 'keep' && featuresHadGaps(fingerprint),
      search: result.search, threshold: result.options.binaryClassificationThreshold,
      thresholdMetric: result.options.thresholdMetric};
    const positiveClass = result.metrics.positiveClass;
    this.results.setAttribute(SHOWN_RUN, runId);
    this.showing++;
    this.publishState(this.problem());
    this.resultsPane.show({engine, table, columns: request.prepared.columns, target: request.target,
      prediction: trainPredictionOf(result), curves: result.curves,
      keptRows: request.keptRows, positiveClass, threshold: result.options.binaryClassificationThreshold,
      blob: result.blob, summary}).catch((e) => reportError(e)).finally(() => {
      this.showing--;
      this.publishState(this.problem());
    });
  }

  /** Empties **Target** and unchecks the **Features** that [table] no longer has; each change is checked as any. */
  private static dropRemovedColumns(table: DG.DataFrame, target: DG.InputBase<DG.Column | null>,
    features: DG.InputBase<DG.Column[]>): void {
    const isGone = (c: DG.Column) => table.col(c.name) === null;
    const chosen = target.value;
    if (chosen !== null && isGone(chosen))
      target.value = null;
    const kept = features.value.filter((c) => !isGone(c));
    if (kept.length < features.value.length)
      features.value = kept;
  }

  /** The inputs of [selection] that a **Runs** row keeps and a restore sets again. */
  private static runInputsOf(selection: TrainingSelection): RunInputs {
    const {engine, target, features, hyperparameters, steps, validation, seed} = selection;
    return {engine, targetName: target.name, featureNames: features.map((c) => c.name), hyperparameters,
      missingValues: selectedMissingValues(selection), steps, twoClasses: twoClasses(target) !== null, validation,
      seed};
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
    return problemMessages(problems);
  }

  private static messageOf(messages: string[] | undefined): string | null {
    return messages?.length ? messages.join(' ') : null;
  }
}
