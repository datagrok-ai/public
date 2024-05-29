import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {category, test, before, expect} from '@datagrok-libraries/utils/src/test';
import {_package, getContainer, getAllModelingEngines} from '../package';
import {readDataframe} from './utils';
import JSZip from 'jszip';

category('chemprop', () => {
  let container: DG.DockerContainer;
  let binBlob: Uint8Array;
  let table: DG.DataFrame;

  before(async () => {
    table = await readDataframe('tests/smiles_test.csv');
    container = await getContainer();
  });

  test('getAllModelingEngines', async () => {
    if (!container) return;
    const modelingEngines = await getAllModelingEngines();
    expect('Chemprop' in DG.toJs(modelingEngines), true);
  });

  test('trainModel', async () => {
    const parameterValues = getParameterValues();
    const tableForPrediction = DG.DataFrame.fromColumns(table.columns.byNames(['canonical_smiles', 'molregno']));
    const modelBlob = await trainModelChemprop(tableForPrediction.toCsv(), 'molregno', parameterValues);

    const zip = new JSZip();
    const archive = await zip.loadAsync(modelBlob);
    const file = archive.file('blob.bin');
    binBlob = await file?.async('uint8array')!;
        
    expect(file !== null, true);
  });

  test('applyModel', async () => {
    const smilesColumn = table.columns.byName('canonical_smiles');
    const column = await applyModelChemprop(binBlob, DG.DataFrame.fromColumns([smilesColumn]).toCsv());
        
    expect(column.length, 29);
  });
});

export async function trainModelChemprop(table: string, predict: string, parameterValues: Record<string, any>): Promise<Uint8Array> {
  const container = await getContainer();
  
  const body = {
    type: 'Chemprop',
    table: table,
    predict: predict,
    parameters: parameterValues
  };

  const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/modeling/train_chemprop', {
    method: 'POST',
    body: JSON.stringify(body),
    headers: { 'Content-Type': 'application/json' }
  });

  if (response.status !== 201)
    throw new Error(`Error training model: ${response.statusText}`);
  return new Uint8Array(await response.arrayBuffer());
}

export async function applyModelChemprop(modelBlob: Uint8Array, table: string): Promise<DG.Column> {
  const container = await getContainer();

  const body = {
    modelBlob: Array.from(modelBlob),
    table: table
  };

  const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/modeling/predict_chemprop', {
    method: 'POST',
    body: JSON.stringify(body),
    headers: { 'Content-Type': 'application/json' }
  });

  if (response.status !== 201)
    throw new Error(`Error applying model: ${response.statusText}`);
  
  const data = await response.json();
  return DG.Column.fromStrings('outcome', data['outcome'].map((v: any) => v?.toString()));
}


function getParameterValues() {
  return {
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
}