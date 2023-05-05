/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';

import {_initEDAAPI} from '../wasm/EDAAPI';
import {computePCA, computePLS} from './EDAtools';
import {renamePCAcolumns, addPLSvisualization} from './EDAui';
import {demoPLS} from './demos';
import {carsDataframe, testDataForBinaryClassification} from './dataGenerators';
import {LINEAR, RBF} from './svm';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//tags: init
export async function init(): Promise<void> {
  await _initEDAAPI();
}

//top-menu: Tools | Data Science | PCA
//name: PCA
//description: Principal component analysis (PCA).
//input: dataframe table
//input: column_list features
//input: int components = 3
//output: dataframe result {action:join(table)}
export async function PCA(table: DG.DataFrame, features: DG.ColumnList, components: number): Promise<DG.DataFrame> {
  return renamePCAcolumns(await computePCA(table, features, components));
}

//top-menu: Tools | Data Science | PLS
//name: PLS
//description: Partial least square regression (PLS).
//input: dataframe table
//input: column_list features
//input: column predict
//input: int components = 3
export async function PLS(table: DG.DataFrame, features: DG.ColumnList, predict: DG.Column, components: number): Promise<void> {
  const plsResults = await computePLS(table, features, predict, components);
  addPLSvisualization(table, features, predict, plsResults);
}

//name: MVA demo
//description: Multivariate analysis (PLS) demo.
//meta.demoPath: Statistical | Multivariate analysis
export async function demoScript(): Promise<any>  {
  const demoScript = new DemoScript('Multivariate analysis', 
    'Provides partial least sqaure regression analysis of the given data.'); 
  
  const cars = carsDataframe();

  const components = 3;
  const predict = cars.columns.byName('price');
  const features = cars.columns.remove('price').remove('model');
  const plsOutput = await computePLS(cars, features, predict, components);  

  const sourceCars = carsDataframe();
  grok.shell.addTableView(sourceCars);
  sourceCars.name = 'Cars';
  let view = grok.shell.getTableView(sourceCars.name);

  await demoScript
    .step('Run', async () => {}, {description: 'Test dataframe is loaded, and multivariate analysis is performed.', delay: 0})
    .step('Study', async () => {addPLSvisualization(sourceCars, features, predict, plsOutput)}, {description: 'Investigate results.', delay: 4000})  
    .step('Try', async () => {
      const params = { items: 10000, features: 100, components: 3};    
      const itemsProp = DG.Property.js('items', DG.TYPE.INT);
      const featuresProp = DG.Property.js('features', DG.TYPE.INT);
      const componentsProp = DG.Property.js('components', DG.TYPE.INT);      
      ui.dialog({title:'Set'})
        .add(ui.input.form(params, [itemsProp, featuresProp, componentsProp]))
        .addButton('Run', async () => await demoPLS(params.items, params.features, params.components))
        .show();      
    }, {description: 'Random walk test dataframe of the given size is generated, and its multivariate analysis is performed.'})    
    .start();
}

//name: Generate linear separable dataset
//description: Generates linear separble dataset for testing binary classificators.
//input: string name = 'Data' {caption: name; category: Dataset}
//input: int samplesCount = 1000 {caption: samples; category: Size}
//input: int featuresCount = 2 {caption: features; category: Size}
//input: double min = -39 {caption: min; category: Range}
//input: double max = 173 {caption: max; category: Range}
//input: double violatorsPercentage = 5 {caption: violators; units: %; category: Dataset}
//output: dataframe df
export async function testDataLinearSeparable(name: string, samplesCount: number, featuresCount: number, 
  min: number, max: number, violatorsPercentage: number): Promise<DG.DataFrame> 
{
  return await testDataForBinaryClassification(LINEAR, [0, 0], name, samplesCount, featuresCount,
    min, max, violatorsPercentage);
}

//name: Generate linear non-separable dataset
//description: Generates linear non-separble dataset for testing binary classificators.
//input: string name = 'Data' {caption: name; category: Dataset}
//input: double sigma = 90  {caption: sigma; category: Hyperparameters} [RBF-kernel paramater]
//input: int samplesCount = 1000 {caption: samples; category: Size}
//input: int featuresCount = 2 {caption: features; category: Size}
//input: double min = -39 {caption: min; category: Range}
//input: double max = 173 {caption: max; category: Range}
//input: double violatorsPercentage = 5 {caption: violators; units: %; category: Dataset}
//output: dataframe df
export async function testDataLinearNonSeparable(name: string, sigma: number, samplesCount: number, 
  featuresCount: number, min: number, max: number, violatorsPercentage: number): Promise<DG.DataFrame> 
{
  return await testDataForBinaryClassification(RBF, [sigma, 0], name, samplesCount, featuresCount,
    min, max, violatorsPercentage);
}
