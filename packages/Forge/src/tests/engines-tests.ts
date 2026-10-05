import * as grok from 'datagrok-api/grok';
import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {defaultHyperparameters, defaultValuesOf, hyperparametersOf, isComplete, rolesOf} from '../engines/engine';
import {EngineRegistry} from '../engines/engine-registry';
import {isApplicable, isInteractive} from '../engines/engine-calls';
import {engineByName, IRIS} from './test-data';

const EDA_ENGINES = ['Linear Regression', 'Softmax', 'PLS Regression', 'XGBoost', 'SVM'];

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
    const without = (...names: string[]) => iris.clone(null, iris.columns.names().filter((n) => !names.includes(n)));
    const species = iris.getCol('Species');
    const petalLength = iris.getCol('Petal.Length');
    const xgboost = engineByName(engines, 'XGBoost');
    const softmax = engineByName(engines, 'Softmax');
    expect(await isApplicable(xgboost, without('Species'), species), true);
    expect(await isInteractive(xgboost, without('Species'), species), true);
    expect(await isApplicable(softmax, without('Species', 'Petal.Length'), petalLength), false);
  });
});
