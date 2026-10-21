import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, expectExceptionAsync, test} from '@datagrok-libraries/test/src/test';
import {applicableEngines} from '../engines/applicable-engines';
import {selectBestEngine} from '../engines/best-engine';
import {defaultHyperparameters, defaultValuesOf, Engine, hasTrainingRows, hyperparametersOf, isComplete,
  isServerEngine, rolesOf} from '../engines/engine';
import {EngineRegistry} from '../engines/engine-registry';
import {isApplicable, isInteractive} from '../engines/engine-calls';
import {SettingDetails} from '../engines/method-list';
import {runListMethods, runTrainModel} from '../training/train-api';
import {columnsOf, deleteRunsOf, demoDemog, engineByName, expectReleased, failsWith, framesSharing, inDiscoveryOrder,
  IRIS, MEASUREMENTS, openDemoTable, valuesOf} from './test-data';

const EDA_ENGINES = ['Linear Regression', 'Softmax', 'PLS Regression', 'XGBoost', 'SVM'];
const BEST_NAMES = ['Chemprop', 'XGBoost', 'PLS Regression', 'Linear Regression'];
const METHOD_COLUMNS = ['Method', 'Package', 'Method type', 'Runs on', 'Keeps training rows', 'Retrains live', 'Roles',
  'Hyperparameters', 'Settings'];

/** An engine with a name only, for the rules that read names. */
function namedEngine(name: string, functions: Engine['functions'] = {}): Engine {
  return {name, namespace: 'Forge', kind: 'function', functions, isLiveUpdate: true};
}

/** A script method's train function, built in the browser (nothing is saved). */
function scriptTrain(language: string, meta: {[key: string]: string} = {}): DG.Script {
  const script = DG.Script.create(`//name: forge-test-train\n//language: ${language}\n` +
    '//input: dataframe df\n//input: string predictColumn\n//output: blob model\n');
  for (const [key, value] of Object.entries(meta))
    script.options[key] = value;
  return script;
}

function floats(name: string, value: (i: number) => number): DG.Column {
  return DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, name, valuesOf(12, value));
}

function strings(name: string, value: (i: number) => string): DG.Column {
  return DG.Column.fromStrings(name, Array.from({length: 12}, (_, i) => value(i)));
}

/** The EDA methods that `listMethods` marks true in [column] for [table] and [target], sorted. */
async function markedMethods(table: DG.DataFrame, target: string | null, column: string): Promise<string[]> {
  const methods = await runListMethods(table, target);
  const names: string[] = methods.getCol('Method').toList();
  return names.filter((name, i) => EDA_ENGINES.includes(name) && methods.get(column, i) === true).sort();
}

category('Engines', () => {
  test('discovers the EDA engines', async () => {
    const engines = EngineRegistry.discover();
    for (const name of EDA_ENGINES) {
      const engine = engineByName(engines, name);
      expect(engine.kind, 'function');
      expect(engine.namespace, 'Eda');
      expect(isComplete(engine), true, `${name} is incomplete`);
    }
  });

  test('groups roles by mlname', async () => {
    const engines = EngineRegistry.discover();
    expect(rolesOf(engineByName(engines, 'PLS Regression')).includes('visualize'), true);
    expect(rolesOf(engineByName(engines, 'XGBoost')).includes('visualize'), false);
    for (const name of EDA_ENGINES) {
      const engine = engineByName(engines, name);
      expect(rolesOf(engine).includes('isInteractive'), true, `${name} has no isInteractive`);
      expect(engine.isLiveUpdate, true, `${name} is not live`);
    }
  });

  test('hyperparameters exclude the table and the target', async () => {
    const xgboost = engineByName(EngineRegistry.discover(), 'XGBoost');
    expectArray(hyperparametersOf(xgboost).map((p) => p.name), ['iterations', 'eta', 'maxDepth', 'lambda', 'alpha']);
  });

  test('default hyperparameters are the train inputs\' initial values', async () => {
    const xgboost = engineByName(EngineRegistry.discover(), 'XGBoost');
    expect(JSON.stringify(defaultHyperparameters(xgboost)),
      JSON.stringify({iterations: 20, eta: 0.3, maxDepth: 6, lambda: 1, alpha: 0}));
  });

  test('a string default loses the quotes the server keeps it in', async () => {
    const kernel = hyperparametersOf(engineByName(EngineRegistry.discover(), 'SVM')).find((p) => p.name === 'kernel');
    if (kernel === undefined)
      throw new Error('SVM has no kernel input');
    expect(/^["']RBF["']$|^RBF$/.test(kernel.initialValue), true, `Initial value ${kernel.initialValue}`);
    expect(defaultValuesOf([kernel])['kernel'], 'RBF');
  });

  test('isApplicable and isInteractive through the contract', async () => {
    const engines = EngineRegistry.discover();
    const iris = await grok.data.files.openTable(IRIS);
    const without = (...names: string[]) => iris.columns.toList().filter((c) => !names.includes(c.name));
    const species = iris.getCol('Species');
    const petalLength = iris.getCol('Petal.Length');
    const xgboost = engineByName(engines, 'XGBoost');
    const softmax = engineByName(engines, 'Softmax');
    expect(await isApplicable(xgboost, without('Species'), species), true);
    expect(await isInteractive(xgboost, without('Species'), species), true);
    expect(await isApplicable(softmax, without('Species', 'Petal.Length'), petalLength), false);
  });

  test('selectBestEngine follows the built-in rule', async () => {
    const all = BEST_NAMES.map((name) => namedEngine(name));
    const withoutChemprop = all.slice(1);
    const pick = (engines: Engine[], features: DG.Column[], target: DG.Column) =>
      selectBestEngine(engines, features, target)?.name;
    const molecule = strings('smiles', (i) => i % 2 === 0 ? 'CCO' : 'c1ccccc1');
    molecule.semType = DG.SEMTYPE.MOLECULE;
    const x = floats('x', (i) => i);
    const label = strings('label', (i) => i % 2 === 0 ? 'a' : 'b');
    const y = floats('y', (i) => 2 * i);
    const six = Array.from({length: 6}, (_, k) => floats(`f${k}`, (i) => i * (k + 1)));

    expect(pick(all, [molecule, x], y), 'Chemprop');
    expect(pick(withoutChemprop, [molecule, x], y), 'XGBoost', 'Without Chemprop the first engine');
    expect(pick(all, [x, strings('kind', (i) => i < 6 ? 'p' : 'q')], label), 'XGBoost');
    expect(pick(all, six, y), 'PLS Regression');
    expect(pick(all, six.slice(0, 2), y), 'Linear Regression');
    expect(pick(all, [strings('kind', (i) => i < 6 ? 'p' : 'q')], label), 'XGBoost');
    expect(pick([], [x], y) === undefined, true, 'An empty list gives a method');
  });

  test('applicableEngines lists the methods for the target, in discovery order', async () => {
    const engines = EngineRegistry.discover();
    const iris = await grok.data.files.openTable(IRIS);
    const species = iris.getCol('Species');
    const petalLength = iris.getCol('Petal.Length');
    const classifiers = (await applicableEngines(engines, columnsOf(iris, MEASUREMENTS), species)).applicable;
    expectArray(classifiers.map((e) => e.name).filter((name) => EDA_ENGINES.includes(name)),
      inDiscoveryOrder(engines, ['XGBoost', 'SVM', 'Softmax']));
    const features = columnsOf(iris, MEASUREMENTS.filter((name) => name !== 'Petal.Length'));
    const [{applicable, failed}, frames] = await framesSharing(features,
      () => applicableEngines(engines, features, petalLength));
    expectArray(applicable.map((e) => e.name).filter((name) => EDA_ENGINES.includes(name)),
      inDiscoveryOrder(engines, ['XGBoost', 'SVM', 'Linear Regression', 'PLS Regression']));
    expect(failed.length, 0, failed.map((f) => f.engine.name).join(', '));
    expectReleased(frames);
  });

  test('applicableEngines leaves out a method whose check throws', async () => {
    const xgboost = engineByName(EngineRegistry.discover(), 'XGBoost');
    const check = xgboost.functions.isApplicable;
    if (check === undefined)
      throw new Error('XGBoost has no isApplicable');
    const throwing: DG.Func = Object.create(check);
    throwing.apply = async () => {
      throw new Error('forge-test: the check failed');
    };
    const broken: Engine = {...xgboost, name: 'forge-test-broken', functions: {...xgboost.functions,
      isApplicable: throwing}};
    const iris = await grok.data.files.openTable(IRIS);
    const {applicable, failed} = await applicableEngines([broken, xgboost], columnsOf(iris, MEASUREMENTS),
      iris.getCol('Species'));
    expectArray(applicable.map((e) => e.name), ['XGBoost']);
    expectArray(failed.map((f) => f.engine.name), ['forge-test-broken']);
    expect(String(failed[0].error).includes('the check failed'), true, String(failed[0].error));
  });

  test('isServerEngine and hasTrainingRows', async () => {
    const engines = EngineRegistry.discover();
    for (const name of EDA_ENGINES) {
      const engine = engineByName(engines, name);
      expect(isServerEngine(engine), false, `${name} runs on the server`);
      expect(hasTrainingRows(engine), name === 'SVM', `${name} training rows`);
    }
    const python = scriptTrain('python');
    expect(python.language, 'python');
    expect(isServerEngine({...namedEngine('forge-test-py', {train: python}), kind: 'script'}), true, 'Python');
    expect(isServerEngine({...namedEngine('forge-test-js', {train: scriptTrain('javascript')}), kind: 'script'}),
      false, 'JavaScript');
    expect(isServerEngine(namedEngine('forge-test-js', {train: scriptTrain('javascript', {mlserver: 'true'})})),
      true, 'meta.mlserver');
    expect(isServerEngine(namedEngine('Chemprop')), true, 'Chemprop');
    expect(hasTrainingRows(namedEngine('forge-test-js', {train: scriptTrain('javascript', {mlhasrows: 'true'})})),
      true, 'meta.mlhasrows');
    expect(hasTrainingRows(namedEngine('forge-test-none')), false, 'No train function');
  });

  test('listMethods lists the complete methods with their settings', async () => {
    const engines = EngineRegistry.discover();
    const methods = await runListMethods(null, null);
    expect(methods.name, 'Methods');
    expectArray(methods.columns.names(), METHOD_COLUMNS);
    const names: string[] = methods.getCol('Method').toList();
    expectArray(names, inDiscoveryOrder(engines, engines.filter(isComplete).map((e) => e.name)));
    for (const name of EDA_ENGINES) {
      const i = names.indexOf(name);
      expect(methods.get('Runs on', i), 'browser', name);
      expect(methods.get('Keeps training rows', i), name === 'SVM', `${name} training rows`);
    }
    const xgboost = engineByName(engines, 'XGBoost');
    const settings: SettingDetails[] = JSON.parse(methods.get('Settings', names.indexOf('XGBoost')));
    expectArray(settings.map((s) => s.name), hyperparametersOf(xgboost).map((p) => p.name));
    const iterations = settings.find((s) => s.name === 'iterations');
    if (iterations === undefined)
      throw new Error('XGBoost has no iterations setting');
    expect(iterations.type, DG.TYPE.INT);
    expect(iterations.default, 20);
    expect(iterations.min, 1);
    expect(iterations.max, 100);
  });

  test('listMethods marks the applicable and the suggested methods of a table', async () => {
    const iris = await grok.data.files.openTable(IRIS);
    expectArray(await markedMethods(iris, 'Species', 'Applicable'), ['SVM', 'Softmax', 'XGBoost']);
    expectArray(await markedMethods(iris, 'Species', 'Suggested'), ['XGBoost']);
    expectArray((await markedMethods(iris, 'Petal.Length', 'Applicable')).filter((name) => name.endsWith('Regression')),
      ['Linear Regression', 'PLS Regression']);
  });

  test('listMethods takes the table\'s last column as the target by default', async () => {
    const iris = await grok.data.files.openTable(IRIS);
    expect(iris.columns.byIndex(iris.columns.length - 1).name, 'Species');
    expectArray(await markedMethods(iris, null, 'Applicable'), ['SVM', 'Softmax', 'XGBoost']);
    expectArray(await markedMethods(iris, null, 'Suggested'), ['XGBoost']);
  });

  test('listMethods refuses a missing target column', async () => {
    const iris = await grok.data.files.openTable(IRIS);
    await expectExceptionAsync(() => runListMethods(iris, 'Specie').then(() => {}),
      failsWith('The table has no column \'Specie\'.'));
    await expectExceptionAsync(() => runListMethods(DG.DataFrame.create(0), null).then(() => {}),
      failsWith('The table has no columns.'));
  });

  test('listMethods refuses a table as trainModel does', async () => {
    const few = DG.DataFrame.fromColumns([DG.Column.fromFloat32Array('x', new Float32Array([1, 2, 3])),
      DG.Column.fromFloat32Array('y', new Float32Array([2, 4, 6]))]);
    const message = 'Training needs at least 10 rows; the table has 3.';
    await expectExceptionAsync(() => runListMethods(few, 'y').then(() => {}), failsWith(message));
    await expectExceptionAsync(() => runTrainModel({table: few, target: 'y', showProgress: false}).then(() => {}),
      failsWith(message));
  });

  test('listMethods suggests the method trainModel trains with', async () => {
    const tables = [await openDemoTable(IRIS, 'methods'), demoDemog(200)];
    try {
      for (const [table, target] of [[tables[0], 'Species'], [tables[1], 'height']] as const) {
        const methods = await runListMethods(table, target);
        const suggested = methods.getCol('Method').toList().filter((_, i) => methods.get('Suggested', i) === true);
        const {method} = await runTrainModel({table, target, showProgress: false});
        expectArray(suggested, [method]);
      }
    } finally {
      await deleteRunsOf(tables);
    }
  }, {timeout: 60000});

  test('listMethods is registered and callable without inputs', async () => {
    const func = DG.Func.find({package: 'Forge', name: 'listMethods'})[0];
    expectArray(func?.inputs.map((p) => p.name) ?? [], ['table', 'target']);
    expect(func?.outputs[0]?.propertyType, DG.TYPE.DATA_FRAME);
    const methods: unknown = await grok.functions.call('Forge:listMethods', {});
    expect(methods instanceof DG.DataFrame && methods.col('Method') !== null, true, `Result: ${methods}`);
  });
});
