import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {category, expectArray, test, before, expect} from '@datagrok-libraries/utils/src/test';
import { _package, getContainer, getAllModelingEngines} from '../package';
import { readDataframe } from './utils';

category('chemprop', () => {
    let container: DG.DockerContainer;
    let modelBlob: Uint8Array;
    let table: DG.DataFrame;
  
    before(async () => {
      table = await readDataframe('tests/smiles_test.csv');
      try {
        container = await getContainer();
      } catch (err: any) {
        _package.logger.error(err);
      }
    });
  
    test('getAllModelingEngines', async () => {
      if (!container) return;
  
      const modelingEngines = await getAllModelingEngines();
      expect('Chemprop' in DG.toJs(modelingEngines), true);
    });

    test('trainModel', async () => {
      const parameterValues = {
        'dataset_type': 'regression',
        'log_frequency': 10,
        'metric': null,
        'multiclass_num_classes': 3,
        'no_cache': false,
        'activation': 'ReLU',
        'atom_messages': false,
        'batch_size': 50,
        'bias': false,
        'depth': 3,
        'dropout': 0,
        'ensemble_size': 1,
        'epochs': 30,
        'ffn_hidden_size': 300,
        'ffn_num_layers': 2,
        'final_lr': 0.0001,
        'hidden_size': 300,
        'init_lr': 0.0001,
        'max_data_size': null,
        'max_lr': 0.001,
        'no_features_scaling': false,
        'num_folds': 1,
        'seed': 0,
        'show_individual_scores': false,
        'split_sizes': [0.8, 0.1, 0.1],
        'split_type': 'random',
        'test': false,
        'undirected': false,
        'use_compound_names': false,
        'warmup_epochs': 2
      };
      modelBlob = await trainModel('111', table.toCsv(), 'molregno', parameterValues);
      expect(modelBlob !== null, true);
    });

    test('applyModel', async () => {
      const smilesColumn = table.columns.byName('canonical_smiles');
      const column = await applyModel('111', modelBlob, DG.DataFrame.fromColumns([smilesColumn]).toCsv());
    });
  });


async function trainModel(id: string, table: string, predict: string, parameterValues: {[_: string]: any}): Promise<Uint8Array> {
  const container = await getContainer();
  const uriParams = new URLSearchParams({
    'id': id,
    'type': 'Chemprop',
    'table': table,
    'predict': predict
  });
  const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id,
    '/modeling/train_chemprop?' + uriParams, {method: 'POST', body: JSON.stringify(parameterValues),
    headers: {'Content-Type': 'application/json'}});
    if (response.status !== 201)
      throw new Error(response.statusText);
    return new Uint8Array(await response.arrayBuffer());
}

export async function applyModel(id: string, modelBlob: Uint8Array, table: string): Promise<DG.Column> {
    const container = await getContainer();
    const uriParams = new URLSearchParams({
      'id': id,
      'type': 'Chemprop',
      'table': table
    });
    const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id,
      '/modeling/predict_chemprop?' + uriParams, {method: 'POST', body: modelBlob,
      headers: {'Content-Type': 'application/octet-stream'}});
    if (response.status !== 201)
      throw new Error(response.statusText);
    const data = await response.json();
    const column = DG.Column.fromStrings('outcome', Array.from(data['outcome'], (v: any, _) => v?.toString()));
    return column;
}