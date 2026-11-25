import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {after, before, category, expect, expectArray, expectFloat, test, waitUntil, waitUntilAsync}
  from '@datagrok-libraries/test/src/test';
import {loadModel} from '../apply/apply-model';
import {exactMapping} from '../apply/column-matching';
import {PANEL_TABLE_MODELS} from '../constants';
import {CHECK_METHOD} from '../engines/check-method';
import {defaultHyperparameters, enabledRules, Engine, hyperparametersOf, isEnabled, isServerEngine,
  visibleHyperparametersOf} from '../engines/engine';
import {errorMessage, ForgeError} from '../forge-error';
import {forgeDb} from '../generated/db';
import {Curves} from '../metrics/metrics';
import {deleteModel} from '../storage/model-store';
import {CandidateRecord, SearchRecord, searchSpaceOf} from '../training/search';
import {trainingEngines, trainPredictionOf} from '../training/train-model';
import {modelLabels} from '../ui/apply-model-dialog';
import {ButtonGate} from '../ui/button-gate';
import {EvaluationView} from '../ui/evaluation-view';
import {FIND_KEY, FindSettings, MethodEntry} from '../ui/find-settings';
import {ForgeApp} from '../ui/forge-app';
import {TIME, TRIALS} from '../ui/limit-inputs';
import {MetricsSummary, MetricsTable, tinyValueText} from '../ui/metrics-table';
import {ModelCards} from '../ui/model-cards';
import {ForgeModelHandler} from '../ui/model-handler';
import {PreparationInputs} from '../ui/preparation-inputs';
import {PanelKey, resultsLayout} from '../ui/results-layout';
import {buildShown, markCurves, refillShown, releaseShown, ResultsData} from '../ui/results-viewers';
import {SHOWN_RUN, TrainView} from '../ui/train-view';
import {AUTO_TUNE, MANUAL, TuneSettings} from '../ui/tune-settings';
import {ValidationInputs} from '../ui/validation-inputs';
import {balloonWith, cancelAtGate, cancelOnceShown, clearTestGalleryStorage, columnsOf, deferred, deleteRunsOf,
  demoDemog, EngineGate, expectReleased, framesSharing, insertModelRow, openIris, openSolubility, runsOf,
  saveIrisModel, slowEngine, stuckEngine, svm, taskBarEntry, twoClassTable, twoSpeciesIris, unbalancedIris,
  useTestGalleryStorage, valuesOf, withEngine, xgboost} from './test-data';

const WAIT_MS = 20000;
const TIMEOUT = 60000;
const STUCK = 'forge-test-stuck';

/** The installed methods and {@link stuckEngine} [STUCK] answering [answer]. */
const withStuck = (answer: Promise<unknown>, calls?: {count: number}): Engine[] =>
  withEngine(stuckEngine(STUCK, answer, calls), trainingEngines().engines.filter((e) => !isServerEngine(e)));

/** The installed methods with XGBoost still retraining live, its train calls waiting for [gate]. */
const withGatedXgboost = (gate: EngineGate): Engine[] =>
  withEngine({...slowEngine(xgboost(), gate), isLiveUpdate: true});

const expectJson = (actual: unknown, expected: unknown) => expect(JSON.stringify(actual), JSON.stringify(expected));

/** The readiness attributes of the results of [view]: `<aria-busy>/<aria-invalid>`. */
const resultsState = (view: DG.ViewBase, results = view.root.querySelector('.forge-train-results')) =>
  `${results?.getAttribute('aria-busy')}/${results?.getAttribute('aria-invalid')}`;
/** The run whose result the results of [view] show (`SHOWN_RUN`), or null. */
const shownRun = (view: DG.ViewBase) =>
  view.root.querySelector('.forge-train-results')?.getAttribute(SHOWN_RUN) ?? null;
const resultsText = (view: TrainView) => view.root.querySelector('.forge-train-results')?.textContent ?? '';
/** What [view] is doing, for the message of a wait that failed. */
const trainState = (view: TrainView) => `training ${view.isTraining}, results ${resultsState(view)}, method ` +
  `${view.methodInput.value}, last training ${view.lastTraining?.engine.name ?? 'none'}`;
/** [view] is checked, trains nothing and shows what it has. */
const isSettled = (view: TrainView) => !view.isTraining && resultsState(view).startsWith('false/');

/** Opens the Train view on [table] (with the [discovered] methods instead of the installed ones when given) and waits
 * for the training it starts by itself. */
async function openTrainView(table: DG.DataFrame, discovered?: Engine[]): Promise<TrainView> {
  const view = await TrainView.create(table, null, discovered);
  grok.shell.addView(view);
  await awaitTrained(view);
  return view;
}

type Training = NonNullable<TrainView['lastTraining']>;

/** Waits for a training of [view] other than [before] that [isWanted] accepts. */
async function awaitTrained(view: TrainView, before?: Training, isWanted: (t: Training) => boolean = () => true):
  Promise<void> {
  await waitUntil(() => view.lastTraining !== undefined && view.lastTraining !== before &&
    isWanted(view.lastTraining) && isSettled(view), () => `No new training: ${trainState(view)}`, TIMEOUT);
}

/** Applies [change] to [view] and waits for the training it starts: a change drops `lastTraining` at once, and the
 * training comes after the check's debounce. */
async function retrain(view: TrainView, change: () => void, isWanted?: (t: Training) => boolean): Promise<void> {
  const before = view.lastTraining;
  change();
  await awaitTrained(view, before, isWanted);
}

/** Closes [view], waits until no training of it is left, then deletes the runs of [tables]. */
async function closeTrainView(view: TrainView, tables: DG.DataFrame[]): Promise<void> {
  view.close();
  await waitUntil(() => !view.isTraining, 'A training outlives the view', TIMEOUT);
  await deleteRunsOf(tables);
}

const isDisabled = (button: HTMLElement) => button.classList.contains('d4-disabled');
const isShown = (element: HTMLElement) => element.style.display !== 'none';
/** Every input of [view]'s form: the data, the method, the preparation, the output, the validation and the method's
 * settings. */
const formInputsOf = (view: TrainView) => [view.tableInput, view.targetInput, view.featuresInput, view.methodInput,
  ...view.preparationInputs.inputs, ...view.outputInputs.inputs, ...view.validationInputs.inputs,
  ...view.hyperparameterInputs.values()];
const formStates = (view: TrainView) => formInputsOf(view).map((input) => `${input.caption} ${input.enabled}`);
/** The enabled inputs of [view]'s form and Settings, then the enabled controls of Search options. */
const enabledInputs = (view: TrainView) => [...[...formInputsOf(view), ...view.tuneSettings.inputs]
  .filter((input) => input.enabled).map((input) => input.caption),
...view.tuneSettings.controls.filter((c) => !c.disabled).map((c) => c.getAttribute('aria-label') ?? c.className)];
/** **Train** is laid out: its row under the inputs is shown. */
const isTrainShown = (view: TrainView) => view.trainButton.offsetParent !== null;

/** Opens the Train view on [table] with SVM as a gated {@link slowEngine} and [target], trains it ([train]) once Train
 * is offered, cancels the training from the task bar at the gate, runs [check] and closes the view. */
async function cancelledSvmTraining(table: DG.DataFrame, target: string, check: (view: TrainView) => Promise<void>,
  train = (view: TrainView) => view.train()): Promise<void> {
  const gate = new EngineGate();
  const view = await TrainView.create(table, 'SVM', withEngine(slowEngine(svm(), gate)));
  view.targetInput.value = table.getCol(target);
  grok.shell.addView(view);
  try {
    await waitUntil(() => isTrainShown(view) && !isDisabled(view.trainButton),
      () => `Train is not offered for SVM: ${trainState(view)}`, WAIT_MS);
    await cancelAtGate(gate, train(view), () => cancelOnceShown('Training SVM model', WAIT_MS));
    await check(view);
  } finally {
    gate.open();
    await closeTrainView(view, [table]);
  }
}

/** The tooltip overlays the platform keeps over the disabled [button]. */
function overlaysOf(button: HTMLElement): HTMLElement[] {
  const rect = button.getBoundingClientRect();
  return Array.from(document.querySelectorAll<HTMLElement>('.d4-tooltip-overlays > span')).filter((o) => {
    const r = o.getBoundingClientRect();
    return Math.abs(r.left - rect.left) < 1 && Math.abs(r.top - rect.top) < 1;
  });
}

/** The value of the hyperparameter input [name] of [view] becomes [value]. */
function setHyperparameter(view: TrainView, name: string, value: number): void {
  const input = view.hyperparameterInputs.get(name);
  if (input === undefined)
    throw new Error(`No ${name} input`);
  input.value = value;
}

/** The platform grid inside [root] whose table has the column [column]. */
function gridWith(root: HTMLElement, column: string): DG.Grid {
  for (const element of Array.from(root.querySelectorAll('.d4-grid'))) {
    const grid: unknown = DG.toJs(DG.Widget.find(element));
    if (grid instanceof DG.Grid && grid.dataFrame.col(column) !== null)
      return grid;
  }
  throw new Error(`No grid with the column ${column}`);
}

/** The texts of the bullet lines of the [lists] under [root]. */
const bulletsOf = (root: HTMLElement, lists = 'ul') =>
  Array.from(root.querySelectorAll(`${lists} > li`), (li) => li.textContent ?? '');

/** Picks [value] in [select], as the user does. */
function pick(select: HTMLSelectElement, value: string): void {
  select.value = value;
  select.dispatchEvent(new Event('change', {bubbles: true}));
}

/** The tests' own storage of the tune settings: the user's stay as they are, and each category's Train views start
 * from the defaults. */
const TEST_TUNE_STORAGE = 'Forge Tune tests';

function useTestTuneStorage(): void {
  TuneSettings.storage = TEST_TUNE_STORAGE;
  clearTestTuneStorage();
}

function clearTestTuneStorage(): void {
  grok.userSettings.put(TEST_TUNE_STORAGE, {});
}

/** The hooks of a category whose tests open a Train view or a Models view, or read or write the settings they keep:
 * the tests' own storages, emptied after the category. */
function testStorageHooks(): void {
  before(async () => {
    useTestTuneStorage();
    useTestGalleryStorage();
  });
  after(async () => {
    clearTestTuneStorage();
    clearTestGalleryStorage();
  });
}

category('UI', () => {
  testStorageHooks();

  test('the core\'s entry points exist; Forge has no menu, app or view function of its own', async () => {
    for (const name of ['forgeTrain', 'forgeApply', 'forgeModels']) {
      const func = DG.Func.find({package: 'Forge', name})[0];
      if (func === undefined)
        throw new Error(`${name} is not registered`);
      expect(!func.topMenu, true, `Unexpected top menu of ${name}: ${func.topMenu}`);
    }
    const app = DG.Func.find({package: 'Forge', name: 'forgeApp'});
    expect(app.length, 1);
    expectArray(app[0].outputs.map((p) => p.propertyType), [DG.TYPE.VIEW]);
    const inputs = DG.Func.find({package: 'Forge', name: 'applyModel'})[0]?.inputs.map((p) => p.name) ?? [];
    expectArray(inputs, ['model', 'table', 'columnNamesMap', 'showProgress']);
    expect(DG.Func.find({package: 'Forge', tags: [DG.FUNC_TYPES.APP]}).length, 0);
    expect(DG.Func.find({package: 'Forge', tags: ['view']}).length, 0);
  });

  test('a change during a training supersedes it without a run row', async () => {
    const iris = await openIris();
    const gate = new EngineGate();
    const view = await TrainView.create(iris, null, withGatedXgboost(gate));
    grok.shell.addView(view);
    try {
      await waitUntil(() => view.isTraining, () => `The view does not train by itself: ${trainState(view)}`, WAIT_MS);
      await retrain(view, () => {
        view.targetInput.value = iris.getCol('Petal.Length');
        view.featuresInput.value = columnsOf(iris, ['Sepal.Length', 'Sepal.Width', 'Petal.Width']);
        gate.open();
      }, (t) => t.result.task === 'regression');
      const recorded = await runsOf([iris]);
      expectArray(recorded.map((r) => `${r.status} ${r.target_name}`), ['completed Petal.Length']);
    } finally {
      gate.open();
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('a cancelled training is recorded as cancelled and leaves nothing to save', async () => {
    const table = twoClassTable(100);
    await cancelledSvmTraining(table, 'label', async (view) => {
      await waitUntil(() => !view.isTraining, 'The cancelled training does not end', WAIT_MS);
      expect(resultsState(view), 'false/true');
      expectArray((await runsOf([table])).map((r) => `${r.status} ${r.engine_name}`), ['cancelled SVM']);
      expect(view.lastTraining === undefined && isDisabled(view.saveButton), true, 'A cancelled training is saveable');
    });
  }, {timeout: TIMEOUT});

  test('Train view drops the old form\'s subscriptions on a Table change', async () => {
    const tables = [grok.shell.addTable(await openIris()), grok.shell.addTable(await openIris())];
    const view = await openTrainView(tables[0]);
    try {
      const oldData = view.groups[0];
      const oldFeatures = view.featuresInput;
      const subCount = view.subs.length;
      for (const table of [tables[1], tables[0]]) {
        const features = view.featuresInput;
        view.tableInput.value = table;
        await waitUntil(() => view.featuresInput !== features, 'The form was not rebuilt', WAIT_MS);
      }
      expect(view.subs.length, subCount);

      const failing = () => 'forge-test-problem';
      for (const [group, input] of [[oldData, oldFeatures], [view.groups[0], view.featuresInput]] as const) {
        group.setExpanded(false);
        input.addValidator(failing);
        input.validate();
      }
      await waitUntil(() => view.groups[0].isExpanded, 'The current Data group does not open on an error', WAIT_MS);
      expect(oldData.isExpanded, false, 'A replaced form still opens its group on an error');
    } finally {
      await closeTrainView(view, tables);
      for (const table of tables)
        grok.shell.closeTable(table);
    }
  }, {timeout: TIMEOUT});

  test('a disabled ButtonGate keeps one overlay', async () => {
    let reason = 'forge-test-reason-0';
    const button = ui.button('forge-test-gate', () => {});
    const gate = new ButtonGate(button, () => reason);
    document.body.append(button);
    try {
      for (let i = 1; i <= 5; i++) {
        gate.update();
        reason = `forge-test-reason-${i}`;
      }
      gate.update();
      // The overlay takes the button's place on the platform's next 100 ms tick; reads 150 ms apart span a later tick.
      let last = -1;
      await waitUntil(() => {
        const count = overlaysOf(button).length;
        const held = count > 0 && count === last;
        last = count;
        return held;
      }, () => `The overlays over the disabled button do not settle: ${last}`, WAIT_MS, 150);
      expect(overlaysOf(button).length, 1);
    } finally {
      button.remove();
    }
  });
});

category('UI: Train', () => {
  testStorageHooks();

  test('Forge Check: checking a molecule feature switches Method to the method that takes it', async () => {
    const table = await openSolubility();
    const view = await openTrainView(table);
    try {
      view.featuresInput.value = [table.getCol('smiles')];
      view.targetInput.value = table.getCol('SOL');
      await waitUntil(() => view.methodInput.value === CHECK_METHOD, `Method is ${view.methodInput.value}`, WAIT_MS);
      await waitUntilAsync(async () => (await runsOf([table])).some((r) => r.engine_name === CHECK_METHOD &&
        r.status === 'completed'), `No completed run of ${CHECK_METHOD}`, TIMEOUT / 2);
      expect(view.methodInput.items.includes('XGBoost'), false, view.methodInput.items.join(', '));
      expect(isShown(view.preparationInputs.group.root), false, 'Preparation is shown');
      expect(isShown(view.outputInputs.group.root), false, 'Output is shown');
      expect(balloonWith('cannot learn') === null, true, 'A balloon named a refused method');
    } finally {
      await closeTrainView(view, [table]);
    }
  }, {timeout: TIMEOUT});

  test('Forge Check: a chosen method that does not take a molecule feature says so', async () => {
    const table = await openSolubility();
    const view = await TrainView.create(table, 'XGBoost');
    grok.shell.addView(view);
    const refusal = `XGBoost cannot learn from Molecule features. Uncheck: smiles. ${CHECK_METHOD} is chosen for now.`;
    try {
      await awaitTrained(view);
      view.featuresInput.value = [table.getCol('smiles'), table.getCol('MolLogP')];
      await waitUntil(() => balloonWith(refusal) !== null, 'No balloon for the refused method', WAIT_MS);
      await waitUntil(() => view.methodInput.value === CHECK_METHOD, `Method is ${view.methodInput.value}`, WAIT_MS);
    } finally {
      await closeTrainView(view, [table]);
    }
  }, {timeout: TIMEOUT});

  test('a training superseded while its task-bar progress is shown writes no run', async () => {
    const iris = await openIris();
    const gate = new EngineGate();
    const view = await TrainView.create(iris, null, withGatedXgboost(gate));
    grok.shell.addView(view);
    const progressShown = () => taskBarEntry('Training XGBoost model') !== undefined;
    try {
      await waitUntil(() => view.isTraining && progressShown(),
        () => `XGBoost does not train with a progress: ${trainState(view)}`, WAIT_MS);
      expect(resultsText(view).includes('Training XGBoost'), false, 'Results names the training');
      await retrain(view, () => {
        setHyperparameter(view, 'maxDepth', 2);
        gate.open();
      });
      const recorded = await runsOf([iris]);
      expectArray(recorded.map((r) => `${r.status} ${r.engine_name} ${r.hyperparameters?.maxDepth}`),
        ['completed XGBoost 2']);
      expect(shownRun(view), recorded[0].id, 'The results show another run than the one recorded');
    } finally {
      gate.open();
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('a cancel unlocks the form; a training completed after a cancelled one clears the invalid mark', async () => {
    const demog = demoDemog(200);
    let before: string[] = [];
    await cancelledSvmTraining(demog, 'sex', async (view) => {
      expect(`${view.runs[0]?.status} ${resultsState(view)}`, 'cancelled false/true',
        'The cancel came after the training ended');
      expectArray(formStates(view), before);
      await view.train();
      await waitUntil(() => resultsState(view).startsWith('false/'), () => `The results stay busy: ${trainState(view)}`,
        WAIT_MS);
      expect(`${view.runs[0]?.status} ${resultsState(view)}`, 'completed false/false');
      expect(shownRun(view), view.lastTraining?.runId ?? null, 'The results show another run than the last training');
    }, (view) => {
      before = formStates(view);
      const training = view.train();
      expectArray(enabledInputs(view), []);
      return training;
    });
  }, {timeout: TIMEOUT});

  test('ValidationInputs: the schemes of the target\'s task, Folds or Test fraction, Split seed with Shuffle',
    async () => {
      const v = new ValidationInputs(() => {});
      expect(v.group.isExpanded, false, 'Validation starts expanded');
      expectArray(v.scheme.items, ['K-fold, stratified', 'K-fold', 'Hold-out']);
      expectArray([v.scheme.value, v.folds.value, v.testFraction.value, v.shuffle.value, v.seed.value],
        ['K-fold, stratified', 5, 0.2, true, 42]);
      expect(isShown(v.folds.root) && !isShown(v.testFraction.root), true, 'Folds is hidden or Test fraction shown');

      v.shuffle.value = false;
      expect(`${v.seed.enabled} ${v.settings().shuffle}`, 'false false', 'Split seed is enabled without Shuffle');
      v.shuffle.value = true;
      expect(v.seed.enabled, true, 'Split seed stays disabled with Shuffle');

      v.scheme.value = 'Hold-out';
      expect(!isShown(v.folds.root) && isShown(v.testFraction.root), true, 'Folds is shown for Hold-out');
      expect(v.settings().scheme, 'holdout');
      v.scheme.value = 'K-fold, stratified';
      v.update(DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'y', [1, 2, 3, 4]));
      expectArray(v.scheme.items, ['K-fold', 'Hold-out']);
      expect(`${v.scheme.value} ${v.settings().stratified}`, 'K-fold false');
      v.update(DG.Column.fromStrings('label', ['a', 'b', 'a', 'b']));
      expectArray(v.scheme.items, ['K-fold, stratified', 'K-fold', 'Hold-out']);
      expect(v.scheme.value, 'K-fold, stratified');

      v.restore({...v.settings(), scheme: 'holdout'}, 7);
      expectArray([v.scheme.value, v.seed.value], ['Hold-out', 7]);
    });

  test('enabledRules of EDA\'s SVM: the settings each kernel and task use', async () => {
    const engine = svm();
    const rules = enabledRules(engine) ?? {};
    const defaults = defaultHyperparameters(engine);
    const states = (kernel: string, task: 'classification' | 'regression', names: string[]) =>
      names.map((name) => `${name} ${isEnabled(rules, name, {...defaults, kernel}, task)}`);
    expect(visibleHyperparametersOf(engine).some((p) => p.name === 'seed'), false, 'SVM has a seed');
    expectArray(states('RBF', 'classification', ['gamma', 'degree', 'coef0', 'epsilon']),
      ['gamma true', 'degree false', 'coef0 false', 'epsilon false']);
    expectArray(states('Polynomial', 'classification', ['degree', 'coef0']), ['degree true', 'coef0 true']);
    expectArray(states('Linear', 'classification', ['gamma', 'degree', 'coef0']),
      ['gamma false', 'degree false', 'coef0 false']);
    expectArray(states('RBF', 'regression', ['epsilon']), ['epsilon true']);
  });

  test('a server method that does not answer does not hold back the training', async () => {
    const iris = await openIris();
    const view = await TrainView.create(iris, null, withStuck(new Promise(() => {})));
    view.checkTimeoutMs = 3000;
    grok.shell.addView(view);
    try {
      await awaitTrained(view);
      expect(view.methodInput.items.includes(STUCK), false, 'Method lists the method still checking');
      expectArray((await runsOf([iris])).map((r) => r.status), ['completed']);
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('a server method\'s late answer joins Method without retraining', async () => {
    const iris = await openIris();
    const answer = deferred<boolean>();
    const view = await TrainView.create(iris, null, withStuck(answer.promise));
    grok.shell.addView(view);
    try {
      await awaitTrained(view);
      const shown = shownRun(view);
      answer.resolve(true);
      await waitUntil(() => view.methodInput.items.includes(STUCK), 'The late method does not join Method', WAIT_MS);
      expect(view.isTraining, false, 'The late answer trains');
      expect(shownRun(view), shown, 'The results show another run');
      expect((await runsOf([iris])).length, 1, 'The late answer recorded a run');
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('a method that does not answer is not asked again', async () => {
    const iris = await openIris();
    const calls = {count: 0};
    const view = await TrainView.create(iris, null, withStuck(new Promise(() => {}), calls));
    view.checkTimeoutMs = 2000;
    grok.shell.addView(view);
    try {
      await awaitTrained(view);
      await waitUntil(() => view.silentMethods.has(STUCK), `${STUCK} is not marked silent`, WAIT_MS);
      await retrain(view,
        () => view.featuresInput.value = columnsOf(iris, ['Sepal.Length', 'Petal.Length', 'Petal.Width']));
      expect(calls.count, 1, 'The silent method was asked again');
      expect(balloonWith('did not answer') === null, true, 'A balloon names the silent method');
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('a chosen server method waits for its answer', async () => {
    const iris = await openIris();
    const answer = deferred<boolean>();
    const view = await TrainView.create(iris, STUCK, withStuck(answer.promise));
    grok.shell.addView(view);
    try {
      await waitUntil(() => resultsState(view) === 'true/false' &&
        resultsText(view).includes(`Waiting for ${STUCK} to check the data...`),
      () => `Results: ${resultsState(view)} ${resultsText(view)}`, WAIT_MS);
      expect(view.methodInput.value, STUCK);
      expect(view.lastTraining === undefined, true, 'The view trained before the method answered');
      expect((await runsOf([iris])).length, 0, 'A run was recorded before the method answered');
      answer.resolve(false);
      await waitUntil(() => balloonWith(`${STUCK} cannot learn from this selection; XGBoost is chosen for now.`) !==
        null, 'No balloon for the refusing method', WAIT_MS);
      await awaitTrained(view);
      expect(view.methodInput.value, 'XGBoost');
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});
});

/** A text column of [rows] values made by [value]. */
const textOf = (rows: number, name: string, value: (i: number) => string) =>
  DG.Column.fromStrings(name, Array.from({length: rows}, (_, i) => value(i)));

category('UI: Preparation', () => {
  test('PreparationInputs: Missing values for gaps, One-hot encoding for text, Skip unique categories for ids',
    async () => {
      const rows = 40;
      const float = (name: string, value: (i: number) => number, nullRows: number[] = []) =>
        DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, name, valuesOf(rows, value, nullRows));
      const height = float('height', (i) => 150 + i, [3]);
      const weight = float('weight', (i) => 60 + i % 9);
      const sex = textOf(rows, 'sex', (i) => i % 2 === 0 ? 'M' : 'F');
      const race = textOf(rows, 'race', (i) => ['A', 'B', 'C'][i % 3]);
      const id = textOf(rows, 'id', (i) => `s${i}`);
      const preparation = new PreparationInputs(() => {}, () => null, []);
      const {oneHot, skipUniqueCategories} = preparation;
      /** The group, Missing values, One-hot encoding and Skip unique categories: shown or not. */
      const shown = () => [isShown(preparation.group.root), isShown(preparation.missingValues.choice.root),
        oneHot.visible, skipUniqueCategories.visible].join(' ');

      preparation.update([height, weight]);
      expect(shown(), 'true true false false', 'Numerical features with a gap');
      preparation.update([weight]);
      expect(shown(), 'false false false false', 'A numerical feature without gaps');
      preparation.update([height, weight, sex, race]);
      expect(shown(), 'true true true false', 'Text features of few categories');
      expect(oneHot.value, true, 'One-hot encoding is not checked for SEX and RACE');
      preparation.update([height, weight, sex, race, id]);
      expect(shown(), 'true true true true', 'An id among the features');
      expectJson(preparation.steps(), {oneHot: true, skipUniqueCategories: true});
      preparation.setKeepsMissingValues(true);
      expect(isShown(preparation.missingValues.choice.root), false, 'Missing values is shown for a method that keeps ' +
        'the gaps');
      preparation.setKeepsMissingValues(false);

      preparation.update([height, weight, sex, race]);
      expect(skipUniqueCategories.visible, false, 'Skip unique categories stays without an id');
      preparation.update([height, weight]);
      expect(`${oneHot.visible} ${oneHot.value}`, 'false false', 'A hidden One-hot encoding is not unchecked');
    });

  test('PreparationInputs: One-hot encoding is checked by default only for few categories and keeps the user\'s ' +
    'choice', async () => {
    const rows = 60;
    const table = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'x', valuesOf(rows, (i) => i)),
      textOf(rows, 'few', (i) => i % 2 === 0 ? 'a' : 'b'), textOf(rows, 'many', (i) => `k${i % 30}`),
      textOf(rows, 'id', (i) => `id${i}`),
    ]);
    const preparation = new PreparationInputs(() => {}, () => null, []);
    const {oneHot, skipUniqueCategories} = preparation;
    const use = (names: string[], isChecked: boolean) => {
      preparation.update(columnsOf(table, names));
      expect(`${oneHot.visible} ${oneHot.value}`, `true ${isChecked}`, `One-hot encoding for ${names.join(', ')}`);
    };
    use(['x', 'few'], true);
    use(['x', 'few', 'many'], false);
    use(['x', 'few', 'id'], true);
    skipUniqueCategories.value = false;
    expect(oneHot.value, false, 'An unskipped id leaves One-hot encoding checked');
    skipUniqueCategories.value = true;
    expect(oneHot.value, true, 'Skipping the id does not check One-hot encoding');

    use(['x', 'many'], false);
    oneHot.value = true;
    use(['x', 'few', 'many'], true);
    preparation.update(columnsOf(table, ['x']));
    expect(`${oneHot.visible} ${oneHot.value}`, 'false false', 'A hidden One-hot encoding is not unchecked');
    use(['x', 'few', 'many'], false);
  });
});

/** The Results rows of [view], in order. */
const rowsOf = (view: TrainView | EvaluationView) =>
  Array.from(view.resultsPane.root.querySelectorAll<HTMLElement>('.forge-results-row'));

const titleOf = (panel: Element) => (panel.firstElementChild?.textContent ?? '').trim();

/** The panel titles of each Results row of [view], in order. */
const rowTitles = (view: TrainView | EvaluationView) =>
  rowsOf(view).map((row) => Array.from(row.querySelectorAll('.forge-results-panel'), titleOf));

/** The viewer of the Results panel titled [title] of [view]. */
function panelViewer(view: TrainView, title: string): DG.Viewer {
  const panel = Array.from(view.resultsPane.root.querySelectorAll<HTMLElement>('.forge-results-panel'))
    .find((p) => titleOf(p) === title);
  const element = panel?.querySelector('.d4-viewer');
  const viewer: unknown = element instanceof HTMLElement ? DG.toJs(DG.Widget.find(element)) : null;
  if (!(viewer instanceof DG.Viewer))
    throw new Error(`No viewer in the panel ${title}`);
  return viewer;
}

/** The row of the metric [label] in a metrics grid. */
const metricRow = (grid: DG.Grid, label: string) => grid.dataFrame.getCol('Metric').toList().indexOf(label);

/** ROC and PR curves over the thresholds 0.75, 0.5, 0.25 and 0.125 (exact in float32); [shift] moves a point of
 * each. */
function curvesOf(shift: number): Curves {
  const values = (list: number[]) => new Float32Array(list);
  return {
    roc: {threshold: values([Infinity, 0.75, 0.5, 0.25, 0.125]), fpr: values([0, 0, 0.25, 0.5, 1]),
      tpr: values([0, 0.5 - shift, 0.75, 1, 1])},
    pr: {threshold: values([0.125, 0.25, 0.5, 0.75, Infinity]), recall: values([1, 1, 0.75, 0.5 - shift, 0]),
      precision: values([0.5, 0.6, 0.75, 1, 1])},
  };
}

/** The results of a model that predicts [target] of [table] (features `a` and `b`) as it is. */
function resultsDataOf(table: DG.DataFrame, target: DG.Column, summary: MetricsSummary, curves?: Curves,
  threshold?: number): ResultsData {
  return {engine: xgboost(), table, columns: columnsOf(table, ['a', 'b']), target,
    prediction: DG.Column.fromStrings('prediction', target.toList().map(String)), curves, keptRows: null, threshold,
    blob: new Uint8Array(0), summary};
}

/** A viewer that throws while it is built, then the Wrong predictions grid: each time the frames go back and no result
 * is shown. */
async function expectFailedBuildReleased(): Promise<void> {
  // One iteration of depth 1 gets rows wrong.
  const unbalanced = await unbalancedIris();
  const viewerFailure = 'forge-test-viewer-failure';
  const gridFailure = 'forge-test-grid-failure';
  const {confusionMatrix, grid} = DG.Viewer;
  let view: TrainView | undefined;
  try {
    DG.Viewer.confusionMatrix = () => {
      throw new ForgeError(viewerFailure);
    };
    const [opened, built] = await framesSharing(unbalanced.columns.toList(), async () => {
      const trainView = await TrainView.create(unbalanced);
      grok.shell.addView(trainView);
      await waitUntil(() => balloonWith(viewerFailure) !== null && !trainView.isTraining,
        () => `The failing viewer was not reported: ${trainState(trainView)}`, TIMEOUT);
      return trainView;
    });
    view = opened;
    expect(view.resultsPane.frames.length, 0, 'A result was shown without its viewers');
    expectReleased(built);

    DG.Viewer.confusionMatrix = confusionMatrix;
    // Only the frame of Wrong predictions, every column of the table and the prediction, gets a grid with both.
    DG.Viewer.grid = (df, options) => {
      if (df.col('Species') !== null && df.col('Predicted Species') !== null)
        throw new ForgeError(gridFailure);
      return grid(df, options);
    };
    const trainView = view;
    const [, rebuilt] = await framesSharing(unbalanced.columns.toList(), async () => {
      setHyperparameter(trainView, 'iterations', 1);
      setHyperparameter(trainView, 'maxDepth', 1);
      await waitUntil(() => balloonWith(gridFailure) !== null && !trainView.isTraining,
        () => `The failing grid was not reported: ${trainState(trainView)}`, TIMEOUT);
    });
    expect(view.resultsPane.frames.length, 0, 'A result was shown without Wrong predictions');
    expectReleased(rebuilt);
  } finally {
    DG.Viewer.confusionMatrix = confusionMatrix;
    DG.Viewer.grid = grid;
    if (view !== undefined)
      await closeTrainView(view, [unbalanced]);
  }
}

category('UI: Results', () => {
  testStorageHooks();

  test('resultsLayout: the template rows, the absent panels dropped, the method\'s rows', async () => {
    const all = new Set<PanelKey>(['predictedVsActual', 'rocCurve', 'prCurve', 'perClass', 'metrics', 'confusions',
      'insights', 'wrongPredictions', 'features', 'correlation', 'classBalance']);
    const binary = resultsLayout('binary', all, 0);
    expectJson(binary.map((r) => r.panels), [['predictedVsActual', 'rocCurve', 'prCurve'],
      ['insights', 'confusions', 'metrics'], ['features', 'correlation', 'classBalance'], ['wrongPredictions']]);
    expectArray(binary.map((r) => r.height), [350, 350, 450, 300]);
    expectJson(binary.map((r) => r.weights), [[1, 1, 1], [1, 1.5, 1.5], [1, 1, 1], [1]]);
    const multiclass = resultsLayout('multiclass', new Set<PanelKey>(['predictedVsActual', 'features', 'correlation',
      'classBalance']), 0);
    expectJson(multiclass.map((r) => r.panels), [['predictedVsActual'], ['insights', 'metrics', 'classBalance'],
      ['features', 'correlation']]);
    expectJson(multiclass.map((r) => r.weights), [[1], [1, 1.5, 1.5], [1, 1]]);
    const regression = resultsLayout('regression', new Set<PanelKey>(['predictedVsActual', 'residuals',
      'distributions', 'statistics']), 4);
    expectJson(regression.map((r) => r.panels), [['predictedVsActual', 'residuals', 'distributions'],
      ['insights', 'metrics', 'statistics'], ['visualize:0', 'visualize:1', 'visualize:2'], ['visualize:3']]);
    expectArray(regression.map((r) => r.height), [350, 350, 350, 350]);
    expectJson(regression.map((r) => r.weights), [[1, 1, 1], [1, 1.5, 1.5], [1, 1, 1], [1]]);
    expectJson(resultsLayout('binary', new Set<PanelKey>(['metrics']), 0).map((r) => r.panels),
      [['insights', 'metrics']]);
  });

  test('resultsLayout: Insights & Tips is always in row 2, Class balance in its template row', async () => {
    const withoutWrong = new Set<PanelKey>(['predictedVsActual', 'confusions', 'insights', 'features', 'correlation',
      'classBalance']);
    expectJson(resultsLayout('multiclass', withoutWrong, 0).map((r) => r.panels), [['predictedVsActual', 'confusions'],
      ['insights', 'metrics', 'classBalance'], ['features', 'correlation']]);
    expectJson(resultsLayout('binary', withoutWrong, 0).map((r) => r.panels), [['predictedVsActual'],
      ['insights', 'confusions', 'metrics'], ['features', 'correlation', 'classBalance']]);
    const full = new Set<PanelKey>([...withoutWrong, 'wrongPredictions', 'perClass']);
    expectJson(resultsLayout('multiclass', full, 0).map((r) => r.panels), [['predictedVsActual', 'confusions',
      'perClass'], ['insights', 'metrics', 'classBalance'], ['features', 'correlation'], ['wrongPredictions']]);
  });

  test('a hovered row is linked through the kept rows', async () => {
    const iris = await openIris();
    iris.getCol('Sepal.Width').set(3, null);
    const view = await openTrainView(iris);
    try {
      await retrain(view, () => view.methodInput.value = 'SVM', (t) => t.engine.name === 'SVM');
      await waitUntil(() => view.resultsPane.frames[0]?.rowCount === 149, 'The results show every row', WAIT_MS);
      const frame = view.resultsPane.frames[0];
      frame.mouseOverRowIdx = 5;
      // The kept rows are 0, 1, 2, 4, 5, 6, ...: the 6th is row 6.
      await waitUntil(() => iris.mouseOverRowIdx === 6, () => `The table hovers row ${iris.mouseOverRowIdx}`, WAIT_MS);
      iris.mouseOverRowIdx = 3;
      await waitUntil(() => frame.mouseOverRowIdx === -1, () => `The frame hovers row ${frame.mouseOverRowIdx}`,
        WAIT_MS);
      expect(iris.mouseOverRowIdx, 3, 'The skipped row came back as none');
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('a retrain keeps the grid, the rows and the viewers, and gives Wrong predictions\' frame back when it goes; a ' +
    'Features change rebuilds the viewers; a hint releases the frames and keeps the grid', async () => {
    const iris = await unbalancedIris();
    const view = await openTrainView(iris);
    try {
      const grid = gridWith(view.resultsPane.root, 'Metric').root;
      const rows = rowsOf(view);
      const titles = rowTitles(view);
      const viewerRoots = () => ['Predicted Species vs Actual', 'Features'].map((t) => panelViewer(view, t).root);
      const roots = viewerRoots();
      const [frame] = view.resultsPane.frames;
      await retrain(view, () => setHyperparameter(view, 'maxDepth', 1));
      expect(gridWith(view.resultsPane.root, 'Metric').root === grid, true, 'A retrain built a new metrics grid');
      const after = rowsOf(view);
      const retrained = rowTitles(view);
      for (let i = 0; i < rows.length; i++) {
        const isSame = JSON.stringify(titles[i]) === JSON.stringify(retrained[i]);
        expect(after[i] === rows[i], isSame, `Row ${i + 1}: ${JSON.stringify(titles[i])}`);
      }
      expect(after[0] === rows[0], true, 'Row 1 was rebuilt');
      expect(viewerRoots().every((root, i) => root === roots[i]), true, 'A retrain rebuilt the viewers');
      expect(view.resultsPane.frames[0] === frame && frame.columns.length > 0, true, 'A retrain replaced the frame');
      const training = view.lastTraining;
      if (training === undefined)
        throw new Error('No training');
      expectArray(frame.getCol('Predicted Species').toList().map(String),
        trainPredictionOf(training.result).toList().map(String));

      // One iteration of depth 1 gets rows wrong; 100 of depth 6 get none of these rows wrong.
      await retrain(view, () => setHyperparameter(view, 'iterations', 1));
      await waitUntil(() => view.resultsPane.frames.length === 2, 'Wrong predictions has no frame of its own', WAIT_MS);
      const mismatches = view.resultsPane.frames[1];
      await retrain(view, () => {
        setHyperparameter(view, 'iterations', 100);
        setHyperparameter(view, 'maxDepth', 6);
      });
      await waitUntil(() => view.resultsPane.frames.length === 1, 'Wrong predictions stays', WAIT_MS);
      expect(view.resultsPane.frames[0] === frame && frame.columns.length > 0, true,
        'The retrain replaced the viewers\' frame');
      expectReleased([mismatches]);

      const shown = [...view.resultsPane.frames];
      await retrain(view, () => view.featuresInput.value = view.featuresInput.value.slice(0, 3));
      await waitUntil(() => view.resultsPane.frames[0] !== frame, 'A Features change kept the frame', WAIT_MS);
      expect(viewerRoots().every((root, i) => root !== roots[i]), true, 'A Features change kept the viewers');
      expectArray(shown.map((f) => f.columns.length), shown.map(() => 0));

      const current = [...view.resultsPane.frames];
      view.resultsPane.showHint('');
      expect(rowsOf(view).length, 0, 'The rows stay with a hint');
      expect(view.resultsPane.frames.length, 0, 'Frames stay with a hint');
      expectArray(current.map((f) => f.columns.length), current.map(() => 0));
      expect(resultsText(view).includes('Performance'), true, 'No Performance header with a hint');
      await retrain(view, () => setHyperparameter(view, 'maxDepth', 6));
      expect(gridWith(view.resultsPane.root, 'Metric').root === grid, true, 'The hint dropped the metrics grid');
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('a viewer or the Wrong predictions grid that fails to build gives its frame back', expectFailedBuildReleased,
    {timeout: TIMEOUT});

  test('buildShown: ROC Curve and PR Curve over frames of their own; markCurves moves only their circles; a refill ' +
    'gives them new frames', async () => {
    const table = twoClassTable(20);
    const target = table.getCol('label');
    const summary: MetricsSummary = {metrics: {train: {auc: 0.8, averagePrecision: 0.7}, validation: {}},
      rowCount: 20, skippedRows: 0};
    const curves = curvesOf(0);
    const data = resultsDataOf(table, target, summary, curves, 0.4);
    const shown = buildShown(data, 'binary');
    const cutoffs = (frame: DG.DataFrame) =>
      frame.getCol('Mark').toList().flatMap((mark, i) => mark === 'cutoff' ? [i] : []);
    try {
      const [roc, pr] = shown.curvePanels;
      expectArray(roc.frame.columns.names(), ['Threshold', 'FPR (1 − Specificity)', 'Recall (Sensitivity)', 'Mark']);
      expectArray(pr.frame.columns.names(), ['Threshold', 'Recall (Sensitivity)', 'Precision', 'Mark']);
      expect(roc.frame.rowCount, 5);
      expect(roc.frame.dart !== shown.frame.dart && roc.viewer.dataFrame.dart === roc.frame.dart, true,
        'ROC Curve is not over a frame of its own');
      expect(roc.frame.meta.formulaLines.items.some((line) => line.title === 'Non-informative baseline'), true,
        'ROC Curve has no baseline');
      expect(pr.frame.meta.formulaLines.items.find((line) => line.title === 'Non-informative baseline')?.formula,
        `\${Precision} = ${curves.pr.precision[0].toFixed(4)}`);
      // The ROC's thresholds decrease (the last at or above 0.4 is circled), the PR's increase (the first).
      expectArray([...cutoffs(roc.frame), ...cutoffs(pr.frame)], [2, 2]);

      markCurves(shown, data);
      expectArray([roc.title.textContent, pr.title.textContent, roc.frame.currentRowIdx, pr.frame.currentRowIdx],
        ['ROC Curve (AUC = 0.800)', 'PR Curve (AP = 0.700)', 2, 2]);
      const [rocFrame, prFrame] = [roc.frame, pr.frame];
      markCurves(shown, {...data, threshold: 0.7});
      expect(roc.frame === rocFrame && pr.frame === prFrame, true, 'A threshold change replaced a curve\'s frame');
      expectArray([...cutoffs(roc.frame), ...cutoffs(pr.frame), roc.frame.currentRowIdx, pr.frame.currentRowIdx],
        [1, 3, 1, 3]);
      markCurves(shown, {...data, threshold: undefined});
      expectArray([...cutoffs(roc.frame), ...cutoffs(pr.frame), roc.frame.currentRowIdx, pr.frame.currentRowIdx],
        [-1, -1]);

      expect(refillShown(shown, resultsDataOf(table, target, summary, curvesOf(0.25)), 'binary'), true,
        'A result of the same data was not refilled');
      expect(roc.frame !== rocFrame && roc.viewer.dataFrame.dart === roc.frame.dart, true,
        'A retrain kept the curve\'s frame');
      expectFloat(roc.frame.get('Recall (Sensitivity)', 1), 0.25, 1e-6);
    } finally {
      releaseShown(shown);
    }
  });

  test('buildShown: Per-class over a frame of its own for three classes, without curves; an axis without values ' +
    'is left out', async () => {
    const rows = 30;
    const table = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'a', valuesOf(rows, (i) => i)),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'b', valuesOf(rows, (i) => (i * 7) % 5)),
      DG.Column.fromStrings('label', Array.from({length: rows}, (_, i) => ['x', 'y', 'z'][i % 3]))]);
    const target = table.getCol('label');
    const summaryOf = (hasAuc: boolean): MetricsSummary => ({metrics: {train: {}, validation: {},
      classes: ['x', 'y', 'z'].map((label, k) => ({label, rows: 10, train: {},
        validation: {precision: 0.9 - k / 10, sensitivity: 0.8, f1: 0.85, ...(hasAuc ? {auc: 0.95} : {})}}))},
    rowCount: rows, skippedRows: 0});
    const perClassOf = (hasAuc: boolean): string[] => {
      const shown = buildShown(resultsDataOf(table, target, summaryOf(hasAuc)), 'multiclass');
      try {
        expect(shown.curvePanels.length, 0, 'Three classes have curves');
        const frame = shown.perClass?.frame;
        if (frame === undefined)
          throw new Error('No Per-class');
        expect(frame.dart !== shown.frame.dart, true, 'Per-class is on the preview frame');
        expectArray(frame.getCol('Class').toList(), ['x', 'y', 'z']);
        expectFloat(frame.get('Precision', 2), 0.7, 1e-6);
        return frame.columns.names();
      } finally {
        releaseShown(shown);
      }
    };
    expectArray(perClassOf(true), ['Class', 'Precision', 'Recall', 'F1 score', 'AUC-ROC']);
    expectArray(perClassOf(false), ['Class', 'Precision', 'Recall', 'F1 score']);
  });

  test('MetricsTable: the Evaluation column and its bullet; a metric only the evaluation has', async () => {
    const evaluated = new MetricsTable({metrics: {train: {accuracy: 0.9, f1: 0.9}, validation: {accuracy: 0.85,
      f1: 0.8}}, rowCount: 150, skippedRows: 0, evaluation: {values: {accuracy: 0.9, f1: 0.9}, tableName: 't',
      rowCount: 150, skippedRows: 0}});
    expectArray(gridWith(evaluated.root, 'Metric').dataFrame.columns.names(),
      ['Metric', 'Train', 'Validation', 'Evaluation']);
    expect(bulletsOf(evaluated.root).includes('Evaluation: t (150 rows)'), true, bulletsOf(evaluated.root).join(' | '));

    const evaluationOnly = gridWith(new MetricsTable({metrics: {train: {accuracy: 0.9}, validation: {accuracy: 0.85}},
      rowCount: 150, skippedRows: 0, evaluation: {values: {accuracy: 0.9, auc: 0.8}, tableName: 't', rowCount: 150,
        skippedRows: 0}}).root, 'Metric');
    const metrics = evaluationOnly.dataFrame;
    expectArray(metrics.getCol('Metric').toList(), ['Accuracy', 'AUC-ROC']);
    const aucRow = metricRow(evaluationOnly, 'AUC-ROC');
    expect(metrics.getCol('Train').isNone(aucRow) && metrics.getCol('Validation').isNone(aucRow), true,
      'AUC-ROC has Train or Validation values');
    expectFloat(metrics.get('Evaluation', aucRow), 0.8, 1e-6);
  });

  test('tinyValueText: a value shown as 0.000 in exponent form', async () => {
    expectArray([0.000032, -0.00004, 0, 0.0004, 0.0012].map(tinyValueText), ['3.2e-05', '-4.0e-05', null, '4.0e-04',
      null]);
  });

  test('EvaluationView: the Evaluation column equals Train on the training data; no column, no application',
    async () => {
      const iris = await openIris();
      const {id} = await saveIrisModel(iris);
      let view: EvaluationView | undefined;
      try {
        const model = await forgeDb.models.get(id);
        const loaded = await loadModel(id);
        const columnCount = iris.columns.length;
        const [opened, frames] = await framesSharing(iris.columns.toList(), () => EvaluationView.open(model, loaded,
          iris, iris.getCol('Species'), exactMapping(loaded.features, iris)));
        view = opened;
        expect(view.name, `Evaluation: ${model.name}`);
        expect(resultsState(view), 'false/null');
        const metrics = gridWith(view.resultsPane.root, 'Metric').dataFrame;
        expectArray(metrics.columns.names(), ['Metric', 'Train', 'Validation', 'Evaluation']);
        for (let i = 0; i < metrics.rowCount; i++)
          expectFloat(metrics.get('Evaluation', i), metrics.get('Train', i), 1e-6);
        const text = view.resultsPane.root.textContent ?? '';
        expect(text.includes(`Evaluation: ${iris.name} (150 rows)`), true, text);
        const titles = rowTitles(view);
        expectJson(titles[0], ['Predicted Species vs Actual', 'Confusions', 'Per-class']);
        expect(titles[1]?.includes('Performance'), true, JSON.stringify(titles));
        expect(titles[1][0], 'Insights & Tips', 'The evaluation hides Insights & Tips');
        // This model gets the overfitting tip, so the metric line stays out.
        expectArray(bulletsOf(view.resultsPane.root, '.forge-insights > ul'), ['The model may be overfitting.']);
        expect(iris.columns.length, columnCount, 'A column was added');
        expect(await forgeDb.applications.query().where('model_id', '=', id).count(), 0);
        view.close();
        view = undefined;
        expectReleased(frames);
      } finally {
        view?.close();
        await deleteModel(id);
      }
    }, {timeout: TIMEOUT});
});

category('UI: Apply dialog', () => {
  test('modelLabels tells apart models with the same name', async () => {
    const t0 = dayjs('2026-01-02T10:00:05');
    const rows = [
      {id: 'a', name: 'N', created_on: t0},
      {id: 'b', name: 'N', created_on: t0.add(1, 'minute')},
      {id: 'c', name: 'N', created_on: t0.add(25, 'second')},
      {id: 'd', name: 'N', created_on: t0.add(25, 'second')},
      {id: 'e', name: 'M', created_on: t0},
    ];
    const labels = [...modelLabels(rows)];
    const labelOf = (id: string) => labels.find(([, row]) => row.id === id)?.[0];
    expect(labels.length, rows.length);
    expect(labelOf('e'), 'M');
    expect(labelOf('b'), 'N (2026-01-02 10:01)');
    expect(labelOf('a'), 'N (2026-01-02 10:00:05)');
    expect(labelOf('c'), 'N (2026-01-02 10:00:30)');
    expect(labelOf('d'), 'N (2026-01-02 10:00:30) #2');
  });
});

category('UI: Catalog', () => {
  test('the autostart registered the model handler', async () => {
    expect(DG.ObjectHandler.list().some((h) => h.name === 'Forge model handler'), true);
    const id = await insertModelRow({name: `forge-test-model-${Date.now()}-handler`});
    try {
      const row = ForgeModelHandler.rowOf(await forgeDb.models.get(id));
      expect(DG.ObjectHandler.forEntity(row)?.name, 'Forge model handler');
    } finally {
      await forgeDb.models.delete(id);
    }
  }, {timeout: TIMEOUT});

  test('Forge:shareModel refuses a missing model without a dialog', async () => {
    const dialogCount = DG.Dialog.getOpenDialogs().length;
    const error = await grok.functions.call('Forge:shareModel', {model: 'forge-test-none'})
      .then(() => null, (e: unknown) => e);
    expect(errorMessage(error).includes('No Forge model \'forge-test-none\' is available to you.'), true,
      String(error));
    expect(DG.Dialog.getOpenDialogs().length, dialogCount, 'A dialog opened');
  }, {timeout: TIMEOUT});
});

category('UI: Table panel', () => {
  testStorageHooks();

  test('tableModelsPanel is a panel of a table', async () => {
    const panel = DG.Func.find({package: 'Forge', name: 'tableModelsPanel'})[0];
    expect(panel?.friendlyName, PANEL_TABLE_MODELS);
    expect(panel?.options['role'], 'panel');
    expect(panel?.options['condition'] === undefined, true, `A condition: ${panel?.options['condition']}`);
    expectArray(panel?.inputs.map((p) => p.propertyType) ?? [], [DG.TYPE.DATA_FRAME]);
  });

  test('the catalog has a description for the AI assistant', async () => {
    const view = await ForgeApp.create();
    try {
      expect(view.aiDescription?.includes('Forge'), true, `aiDescription: ${view.aiDescription}`);
    } finally {
      view.close();
    }
  }, {timeout: TIMEOUT});
});

/** The **Runs** grid of [view]; fails without one. */
function runsGrid(view: TrainView): DG.Grid {
  const grid = view.runsPanel.grid;
  if (grid === null)
    throw new Error('No Runs grid');
  return grid;
}

const runsColumn = (view: TrainView, name: string) => runsGrid(view).dataFrame.getCol(name).toList();

category('UI: Runs', () => {
  testStorageHooks();

  test('restoreRun of a run whose inputs differ sets them and trains; of the form\'s own inputs sets nothing',
    async () => {
      const iris = await openIris();
      const view = await openTrainView(iris);
      try {
        await retrain(view, () => setHyperparameter(view, 'iterations', 30));
        await retrain(view, () => {
          view.restoreRun(view.runs[1]);
          expect(view.hyperparameterInputs.get('iterations')?.value, 20, 'The run\'s inputs are not set');
        });
        expectArray(runsColumn(view, '#'), [3, 2, 1]);
        expect(runsGrid(view).dataFrame.currentRowIdx, 0);

        // A run of other inputs drops the result at once: it stays only when nothing is set.
        const third = view.lastTraining;
        view.restoreRun(view.runs[0]);
        expect(view.lastTraining === third && !view.isTraining, true, 'The form\'s own inputs were set again');

        view.featuresInput.value = [];
        await waitUntil(() => resultsState(view) === 'false/true',
          () => `No problem without features: ${trainState(view)}`, WAIT_MS);
        expect(view.runs.length, 3, 'An empty Features trained');
        view.restoreRun(view.runs[0]);
        expect(view.featuresInput.value.length, 4, 'The run\'s features are not set again');
        await awaitTrained(view);
        await waitUntil(() => view.runs.length === 4, 'A run of other inputs does not train', WAIT_MS);
      } finally {
        await closeTrainView(view, [iris]);
      }
    }, {timeout: TIMEOUT});

  test('restoreRun of the current run after a threshold move restores the threshold without training', async () => {
    const iris = await twoSpeciesIris();
    const view = await openTrainView(iris);
    const {predictProbability, threshold} = view.outputInputs;
    const hasProbabilityStep = () =>
      view.lastTraining?.result.options.postprocessingInfo.includes('class-probabilities');
    try {
      predictProbability.value = true;
      await waitUntil(() => hasProbabilityStep() === true && isSettled(view),
        () => `Predict probability does not re-cut: ${trainState(view)}`, WAIT_MS);
      const [count, shown] = [view.runs.length, shownRun(view)];
      threshold.value = 0.7;
      await waitUntil(() => view.lastTraining?.result.options.binaryClassificationThreshold === 0.7,
        'The threshold does not re-cut', WAIT_MS);
      expect(shown !== null && shownRun(view) === shown, true, `A re-cut shows ${shownRun(view)}, not ${shown}`);
      view.restoreRun(view.runs[0]);
      await waitUntil(() => hasProbabilityStep() === false && isSettled(view),
        () => `The run's Predict probability does not re-cut: ${trainState(view)}`, WAIT_MS);
      expect(predictProbability.value, view.runs[0].steps.predictProbability);
      expect(`${view.runs[0].steps.predictProbability} ${threshold.value}`, 'false 0.7',
        'Unchecking by a run without Predict probability did not keep the hidden threshold');
      expect(!view.isTraining && view.runs.length === count && shownRun(view) === shown, true, 'The restore trained');
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('restoreRun of a run whose columns the table lost sets nothing', async () => {
    const iris = await openIris();
    const view = await openTrainView(iris);
    try {
      const first = view.runs[0];
      await retrain(view, () => iris.columns.remove('Petal.Width'));
      const [count, training] = [view.runs.length, view.lastTraining];
      const features = view.featuresInput.value.map((c) => c.name);
      view.restoreRun(first);
      expect(view.lastTraining === training && !view.isTraining, true, 'A run with a missing column set the form');
      await waitUntil(() => balloonWith('The table no longer has the columns of this run') !== null,
        () => `No balloon for the missing columns: ${trainState(view)}`, WAIT_MS);
      const balloon = balloonWith('The table no longer has the columns of this run');
      expect(balloon?.textContent?.includes('The table no longer has the columns of this run: Petal.Width.'), true,
        balloon?.textContent ?? '');
      expect(view.runs.length, count);
      expectArray(view.featuresInput.value.map((c) => c.name), features);
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});
});

const searchRecord = (holdOut?: {rows: number; fraction: number}): SearchRecord => ({mode: 'tune',
  metric: 'accuracy', budgetSeconds: 60, seed: 1, trials: 12, completed: 10, pruned: 1, failed: 1, elapsedMs: 1000,
  stopped: 'budget', ...(holdOut === undefined ? {} : {holdOut}), leaders: [{method: 'XGBoost',
    hyperparameters: {iterations: 20}, score: 0.8, se: 0.01, foldScores: [0.8, 0.8, 0.8, 0.8, 0.8], trials: 12,
    chosen: true}]});

category('UI: Search', () => {
  testStorageHooks();

  test('TuneSettings lists the metrics of the task, Log loss and Brier score last', async () => {
    const select: HTMLSelectElement = ui.element('select');
    const labels = () => Array.from(select.options, (o) => o.textContent ?? '');
    expect(TuneSettings.listMetrics(select, undefined, 'classification', undefined), 'classification');
    expectArray(labels().slice(-2), ['Log loss', 'Brier score']);
    expect(select.value, 'balancedAccuracy');
    TuneSettings.listMetrics(select, 'classification', 'regression', 'mae');
    expectArray(labels(), ['R2', 'RMSE', 'MAE']);
    expect(select.value, 'mae');
  });

  test('TuneSettings.load: a stale or malformed item falls back to its default; old seconds become whole minutes',
    async () => {
      const engine = xgboost();
      const space = searchSpaceOf(engine);
      const properties = hyperparametersOf(engine);
      const loaded = (value: unknown) => {
        grok.userSettings.add(TuneSettings.storage, 'XGBoost', JSON.stringify(value));
        return TuneSettings.load(engine, space, properties);
      };
      try {
        const kept = loaded({mode: AUTO_TUNE, expanded: true, pinned: ['maxDepth'],
          ranges: {iterations: {low: 15, high: 100}}, metric: 'f1', limit: TRIALS, minutes: 5, trials: 20});
        expectJson([kept.mode, kept.expanded, kept.pinned, kept.ranges, kept.metric, kept.limit, kept.minutes,
          kept.trials], [AUTO_TUNE, true, ['maxDepth'], {iterations: {low: 15, high: 100}}, 'f1', TRIALS, 5, 20]);
        const stale = loaded({mode: AUTO_TUNE, pinned: ['speed', 'eta'], ranges: {iterations: {low: 5, high: 500},
          eta: {low: 0.1, high: 0.2}}, metric: 'r2', trials: 0, seconds: 90});
        expectJson([stale.mode, stale.expanded, stale.pinned, stale.ranges, stale.limit, stale.minutes, stale.trials],
          [AUTO_TUNE, false, ['eta'], {eta: {low: 0.1, high: 0.2}}, TIME, 2, 50]);
        const limits = (value: unknown) => {
          const {limit, minutes, trials} = loaded(value);
          return `${limit} ${minutes} ${trials}`;
        };
        expectArray([limits({seconds: 30}), limits({seconds: 900, trials: 70}), limits({seconds: 7200}),
          limits({limit: TRIALS, minutes: 61, trials: 10001}), limits({limit: TRIALS, minutes: 3, trials: 9}),
          limits({limit: TIME, minutes: 4}), limits({mode: AUTO_TUNE}), limits({limit: 'Steps'})],
        [`${TIME} 1 50`, `${TIME} 15 70`, `${TIME} 60 50`, `${TRIALS} 1 50`, `${TRIALS} 3 9`, `${TIME} 4 50`,
          `${TRIALS} 1 50`, `${TRIALS} 1 50`]);
        grok.userSettings.add(TuneSettings.storage, 'XGBoost', '{not json');
        expectJson(TuneSettings.load(engine, space, properties), {mode: MANUAL, expanded: false, pinned: [],
          ranges: {}, limit: TRIALS, minutes: 1, trials: 50});
      } finally {
        grok.userSettings.delete(TuneSettings.storage, 'XGBoost');
      }
    }, {timeout: TIMEOUT});

  test('MetricsTable after a search: the Test column and the search\'s bullets', async () => {
    const summary: MetricsSummary = {metrics: {train: {accuracy: 0.9}, validation: {accuracy: 0.8},
      test: {accuracy: 0.85}}, rowCount: 150, skippedRows: 0, search: searchRecord({rows: 1170, fraction: 0.2})};
    const table = new MetricsTable(summary);
    expectArray(gridWith(table.root, 'Metric').dataFrame.columns.names(), ['Metric', 'Train', 'Validation', 'Test']);
    const lines = bulletsOf(table.root);
    expect(lines.includes('Search: Tune by Accuracy, limit 1 min: 12 trials of 1 method, chosen by the ' +
      'one-standard-error rule'), true, lines.join(' | '));
    expect(lines.includes('Test: 1170 rows held out before the search'), true, lines.join(' | '));
    const textOf = (search: SearchRecord) => new MetricsTable({...summary, search}).root.textContent ?? '';
    const halving = textOf({...searchRecord(), mode: 'best', budgetSeconds: undefined, maxTrials: 100,
      strategy: 'halving'});
    expect(halving.includes('Search: Find best model by Accuracy, limit 100 trials, Focus on leaders: 1 method, ' +
      'chosen by the highest Accuracy − SE'), true, halving);
    const continued = textOf({...searchRecord(), maxTrials: 50});
    expect(continued.includes('Search: Tune by Accuracy, limit 1 min + 50 trials: 12 trials of 1 method'), true,
      continued);
    const optimistic = textOf(searchRecord());
    expect(optimistic.includes('Validation after a search is optimistic'), true, optimistic);
  });
});

/** Each row of [table] as its cells' texts. */
const tableRows = (table: HTMLElement) =>
  Array.from(table.querySelectorAll('tr'), (row) => Array.from(row.cells, (cell) => cell.textContent ?? ''));

category('UI: Cards', () => {
  test('bestMarks: 3 decimals, ties, all when all tie, not one card; kindLine, settingsOf, settingsBlock', async () => {
    expectArray(ModelCards.bestMarks([0.9801, 0.9804, 0.97], 1), [true, true, false]);
    expectArray(ModelCards.bestMarks([0.03, 0.02, 0.0204], -1), [false, true, true]);
    expectArray(ModelCards.bestMarks([0.5], 1), [false]);
    expectArray(ModelCards.bestMarks([0.5, 0.5], 1), [true, true]);
    expectArray(ModelCards.bestMarks([0.5001, 0.4999], 1), [true, true]);
    expectArray(ModelCards.bestMarks([undefined, undefined], 1), [false, false]);
    expectArray(ModelCards.bestMarks([0.5, undefined], 1), [true, false]);

    const engine = xgboost();
    const defaults = defaultHyperparameters(engine);
    const candidate = (method: string): CandidateRecord => ({method, trials: 4, outcome: 'complete'});
    const find = (methods: string[]): SearchRecord =>
      ({...searchRecord(), mode: 'best', methods: methods.map(candidate)});
    expectArray([ModelCards.kindLine({settingsStatus: 'Default'}), ModelCards.kindLine({settingsStatus: 'Customized'}),
      ModelCards.kindLine({settingsStatus: 'Auto-tuned', search: searchRecord()}),
      ModelCards.kindLine({settingsStatus: 'Auto-tuned', search: find(['XGBoost'])}),
      ModelCards.kindLine({settingsStatus: 'Auto-tuned', search: find(['XGBoost', 'SVM'])})],
    ['Default settings', 'Customized settings', 'Auto-tuned settings', 'Auto-tuned settings', 'Auto-tuned model']);

    const properties = visibleHyperparametersOf(engine);
    const entries = (values: Record<string, unknown>) => properties
      .map((p) => `${DG.StringUtils.camelCaseToSentence(p.caption || p.name)} ${String(values[p.name])}`);
    expect(properties.some((p) => p.name === 'iterations'), true, 'XGBoost has no iterations');
    const shown = (values: typeof defaults) => ModelCards.settingsOf({engine, hyperparameters: values})
      .map((s) => `${s.caption} ${s.value}`);
    expectArray(shown(defaults), entries(defaults));
    expectArray(shown({...defaults, iterations: 30}), entries({...defaults, iterations: 30}));
    expectArray(ModelCards.settingsOf({engine, hyperparameters: defaults}).map((s) => s.description),
      properties.map((p) => p.description ?? ''));

    expectArray(ModelCards.settingsBlock('Default settings', []), []);
    const settings = (n: number) => Array.from({length: n}, (_, i) =>
      ({caption: `Setting ${i}`, value: `${i}`, description: i === 0 ? '' : `What ${i} does`}));
    /** The text of [block]'s kind line, then its one table's rows, each row's cells joined by spaces. */
    const blockOf = (block: HTMLElement[]) => [block[0].textContent ?? '',
      ...tableRows(block[1].querySelector('table')!).map((r) => r.join(' '))];
    expect(ModelCards.settingsBlock('Auto-tuned settings', settings(5))[1].querySelectorAll('table').length, 1);
    expectArray(blockOf(ModelCards.settingsBlock('Auto-tuned model', settings(1))),
      ['Auto-tuned model', 'Setting 0 0']);
    expectArray(blockOf(ModelCards.settingsBlock('Customized settings', settings(3))),
      ['Customized settings', 'Setting 0 0 Setting 1 1', 'Setting 2 2']);
    expectArray(blockOf(ModelCards.settingsBlock('Auto-tuned model', settings(4))),
      ['Auto-tuned model', 'Setting 0 0 Setting 1 1', 'Setting 2 2 Setting 3 3']);
    expectArray(blockOf(ModelCards.settingsBlock('Auto-tuned settings', settings(5))),
      ['Auto-tuned settings', 'Setting 0 0 Setting 1 1', 'Setting 2 2 Setting 3 3', 'Setting 4 4']);
  });
});

category('UI: Find best', () => {
  testStorageHooks();

  test('FindSettings.load: Strategy, the limit and the method toggles; an old stored time becomes whole minutes',
    async () => {
      const loaded = (value: unknown) => {
        grok.userSettings.add(TuneSettings.storage, FIND_KEY, JSON.stringify(value));
        return FindSettings.load();
      };
      const limits = (value: unknown) => {
        const {strategy, limit, minutes, trials} = loaded(value);
        return `${strategy} ${limit} ${minutes} ${trials}`;
      };
      try {
        expectArray([limits({seconds: 30, trials: 70}), limits({strategy: 'halving', limit: TRIALS, minutes: 7,
          trials: 40}), limits({strategy: 'fast', seconds: 3600.5})],
        [`equal ${TIME} 1 70`, `halving ${TRIALS} 7 40`, `equal ${TIME} 60 100`]);
        expectJson(loaded({methods: {XGBoost: false, SVM: 'yes'}}).methods, {XGBoost: false});
      } finally {
        grok.userSettings.delete(TuneSettings.storage, FIND_KEY);
      }
    }, {timeout: TIMEOUT});

  test('FindSettings: Candidates checks the quick methods at a quick limit, every one above it, and lists only the ' +
    'methods the server and probability rules allow', async () => {
    const entry = (name: string, isLive: boolean, probabilities: MethodEntry['probabilities'], isServer = false) =>
      ({name, isLive, isServer, probabilities});
    const noProbabilities = 'forge-test-no-probabilities';
    const find = new FindSettings(ui.button('Find', () => {}), ui.label('Method'), {onChanged: () => {}});
    try {
      find.setTask('classification');
      find.setMethods([entry('XGBoost', true, 'always'), entry('SVM', true, 'input'),
        entry('RealMLP', false, 'always'), entry(noProbabilities, true, 'none'),
        entry('forge-test-server', true, 'always', true)]);
      expectArray(find.candidates.items, ['XGBoost', 'SVM', 'RealMLP', noProbabilities]);
      expectArray(find.checkedMethods, ['XGBoost', noProbabilities]);
      find.limit.kind.value = TIME;
      find.limit.minutes.value = 5;
      expectArray(find.checkedMethods, ['XGBoost', 'SVM', 'RealMLP', noProbabilities]);
      find.limit.minutes.value = 1;
      expectArray(find.checkedMethods, ['XGBoost', noProbabilities]);

      find.setPredictProbability(true);
      expectArray(find.candidates.items, ['XGBoost', 'SVM', 'RealMLP']);
      find.setPredictProbability(false);
      pick(find.metric, 'auc');
      expectArray(find.candidates.items, ['XGBoost', 'SVM', 'RealMLP']);

      find.candidates.value = [];
      expect(find.problem(), 'Check at least one method');
    } finally {
      for (const sub of find.subs)
        sub.unsubscribe();
    }
  }, {timeout: TIMEOUT});
});
