import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {applicableEngines} from '../engines/applicable-engines';
import {selectBestEngine} from '../engines/best-engine';
import {defaultHyperparameters, defaultValuesOf, Engine, hasTrainingRows, hyperparametersOf, isComplete,
  isServerEngine, rolesOf} from '../engines/engine';
import {EngineRegistry} from '../engines/engine-registry';
import {isApplicable, isInteractive} from '../engines/engine-calls';
import {columnsOf, engineByName, expectReleased, framesSharing, inDiscoveryOrder, IRIS, MEASUREMENTS, valuesOf}
  from './test-data';

const EDA_ENGINES = ['Linear Regression', 'Softmax', 'PLS Regression', 'XGBoost', 'SVM'];
const BEST_NAMES = ['Chemprop', 'XGBoost', 'PLS Regression', 'Linear Regression'];

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
});
