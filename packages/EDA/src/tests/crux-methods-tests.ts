import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, expect, expectExceptionAsync, test} from '@datagrok-libraries/test/src/test';

import {isCatBoostApplicable, isCatBoostInteractive, isDecisionTreeApplicable, isDecisionTreeInteractive,
  isExtraTreesApplicable, isExtraTreesInteractive, isHistGradientBoostingApplicable,
  isHistGradientBoostingInteractive, isKnnApplicable, isKnnInteractive, isLdaApplicable, isLdaInteractive,
  isMlpApplicable, isMlpInteractive, isNaiveBayesApplicable, isNaiveBayesInteractive, isRandomForestApplicable,
  isRandomForestInteractive, isRealMlpApplicable, isRealMlpInteractive, isTabmApplicable, isTabmInteractive,
  trainCatBoost, trainHistGradientBoosting, trainKnn, trainLda, trainMlp, trainRandomForest, trainRealMlp,
  trainTabm} from '../crux/methods';
import {readCruxManifest} from '../crux/crux-model';
import {Fit, FITS, loadCase} from './baseline';
import {accuracy, applyTwice, asTable, classificationDataset, copyCol, demogFrame, expectGapsKept, expectRejected,
  GAP_ROWS, irisMeasurements, numericFeatures, predictColumn, roundTrip, sameValues, seededRegression,
  splitTarget} from './utils';

const IRIS = 'System:DemoFiles/iris.csv';

const TIMEOUT = 60000;
const BENCH_TIMEOUT = 600000;
const MAX_BLOB_BYTES = 20 * 1024 * 1024;
const MAX_INTERACTIVE_FIT_MS = 1000;

type Predicate = (features: DG.ColumnList, target: DG.Column) => boolean;
type Shape = {rows: number, features: number, target: 'float' | 'string'};

interface MethodCase {
  name: string;
  fit: Fit;
  regression: boolean;
  classification: boolean;
  nativeNans: boolean;
  standard: boolean;
  applicable: Predicate;
  interactive: Predicate;
  params: {regression?: object, classification?: object};
  floors: {species: number, petal: number};
  seeded: boolean;
  shapes: {applicable: [Shape, boolean][], interactive: [Shape, boolean][]};
  tests?: [string, () => Promise<void>][];
  boolRows?: number;
}

const floatShape = (rows: number, features: number): Shape => ({rows, features, target: 'float'});
const stringShape = (rows: number, features: number): Shape => ({rows, features, target: 'string'});

let irisFile: Promise<DG.DataFrame> | null = null;

/** iris read once per bundle; every test gets its own copy */
async function readIris(): Promise<DG.DataFrame> {
  irisFile ??= grok.dapi.files.readCsv(IRIS);
  return (await irisFile).clone();
}

async function irisSpecies(): Promise<{iris: DG.DataFrame, df: DG.DataFrame, species: DG.Column}> {
  const iris = await readIris();
  return {iris, df: irisMeasurements(iris), species: iris.col('Species')!};
}

const paramsOf = async (bytes: Promise<Uint8Array>) => readCruxManifest(await bytes).params;

const forestParams = (model: string) => ({
  classification: {n_estimators: 100, max_depth: 10, min_samples_leaf: 1, max_features: 'sqrt', random_state: 0,
    model},
  regression: {n_estimators: 100, max_depth: 10, min_samples_leaf: 1, max_features: 'all', random_state: 0, model},
});

const both = (params: object) => ({classification: params, regression: params});

const FOREST_INTERACTIVE: [Shape, boolean][] = [[floatShape(2000, 10), true], [floatShape(2001, 10), false]];

const BOOSTING_APPLICABLE: [Shape, boolean][] = [[floatShape(100000, 100), true], [floatShape(100001, 100), false]];

const CLASSIFIER_APPLICABLE: [Shape, boolean][] = [[stringShape(100000, 100), true], [stringShape(100001, 100), false]];

const TREES = {regression: true, classification: true, nativeNans: true, standard: false, seeded: true,
  floors: {species: 0.9, petal: 0.8}};

const NON_NATIVE = {...TREES, nativeNans: false, seeded: false};

const CLASSIFIER_ONLY = {regression: false, floors: {species: 0.9, petal: 0}};

const METHODS: MethodCase[] = [
  {name: 'Random Forest', fit: FITS.randomForest, ...TREES, applicable: isRandomForestApplicable,
    interactive: isRandomForestInteractive, params: forestParams('random_forest'),
    shapes: {applicable: [[floatShape(100000, 10), true], [floatShape(100001, 10), false]],
      interactive: FOREST_INTERACTIVE},
    tests: [['Features per split', async () => {
      const {iris, df, species} = await irisSpecies();
      const petal = iris.col('Petal.Length')!;
      expect((await paramsOf(trainRandomForest(df, species, 10, 10, 1, 'auto', 0))).max_features, 'sqrt');
      expect((await paramsOf(trainRandomForest(df, petal, 10, 10, 1, 'auto', 0))).max_features, 'all');
      expect((await paramsOf(trainRandomForest(df, petal, 10, 10, 1, 'log2', 0))).max_features, 'log2');
    }]]},
  {name: 'Extra Trees', fit: FITS.extraTrees, ...TREES, applicable: isExtraTreesApplicable,
    interactive: isExtraTreesInteractive, params: forestParams('extra_trees'),
    shapes: {applicable: [[floatShape(200000, 10), true], [floatShape(200001, 10), false]],
      interactive: FOREST_INTERACTIVE}},
  {name: 'Hist Gradient Boosting', fit: FITS.histGradientBoosting, ...TREES,
    applicable: isHistGradientBoostingApplicable, interactive: isHistGradientBoostingInteractive,
    params: both({n_estimators: 100, learning_rate: 0.1, max_leaf_nodes: 31, min_samples_leaf: 20,
      l2_regularization: 0, early_stopping: 'off', random_state: 0, model: 'hist_gradient_boosting'}),
    shapes: {applicable: BOOSTING_APPLICABLE,
      interactive: [[floatShape(20000, 10), true], [floatShape(20001, 1), false], [floatShape(10000, 21), false]]},
    tests: [['Max depth 0 means none', async () => {
      const {df, species} = await irisSpecies();
      expect('max_depth' in await paramsOf(trainHistGradientBoosting(df, species, 10, 0.1, 31, 0, 20, 0, 0)), false,
        'max depth 0 must not be sent');
      expect((await paramsOf(trainHistGradientBoosting(df, species, 10, 0.1, 31, 5, 20, 0, 0))).max_depth, 5);
    }]]},
  {name: 'CatBoost', fit: FITS.catBoost, ...TREES, applicable: isCatBoostApplicable,
    interactive: isCatBoostInteractive,
    params: both({iterations: 100, learning_rate: 0.1, depth: 6, l2_leaf_reg: 3, random_state: 0, model: 'catboost'}),
    shapes: {applicable: BOOSTING_APPLICABLE,
      interactive: [[floatShape(10000, 10), true], [floatShape(20001, 1), false], [floatShape(5000, 21), false]]},
    tests: [['Validation', async () => {
      const {df, species} = await irisSpecies();
      await expectRejected([
        ['iterations 0', () => trainCatBoost(df, species, 0, 0.1, 6, 3, 0)],
        ['depth 17', () => trainCatBoost(df, species, 10, 0.1, 17, 3, 0)],
        ['learning rate 0', () => trainCatBoost(df, species, 10, 0, 6, 3, 0)],
        ['learning rate 1.5', () => trainCatBoost(df, species, 10, 1.5, 6, 3, 0)],
        ['L2 -1', () => trainCatBoost(df, species, 10, 0.1, 6, -1, 0)],
      ]);
    }]]},
  {name: 'Decision Tree', fit: FITS.decisionTree, ...TREES, applicable: isDecisionTreeApplicable,
    interactive: isDecisionTreeInteractive,
    params: both({max_depth: 10, min_samples_leaf: 1, random_state: 0, model: 'decision_tree'}),
    shapes: {applicable: BOOSTING_APPLICABLE,
      interactive: [[floatShape(20000, 10), true], [floatShape(20001, 10), false]]}},
  {name: 'k-Nearest Neighbors', fit: FITS.knn, ...NON_NATIVE, standard: true, applicable: isKnnApplicable,
    interactive: isKnnInteractive, params: both({n_neighbors: 5, weights: 'uniform', model: 'knn'}),
    shapes: {
      applicable: [[floatShape(50000, 50), true], [floatShape(50001, 1), false], [floatShape(25000, 101), false]],
      interactive: [[floatShape(10000, 10), true], [floatShape(10001, 1), false], [floatShape(5000, 21), false]],
    },
    tests: [['Neighbors above the rows', async () => {
      const {iris, df, species} = await irisSpecies();
      await expectExceptionAsync(async () => {
        await trainKnn(df, species, 200, 'Uniform');
      }, (e) => e.message.includes('not above the training rows (150)'));
      const part = iris.clone(DG.BitSet.create(iris.rowCount, (i) => i < 120));
      await expectExceptionAsync(async () => {
        await trainKnn(irisMeasurements(part), part.col('Species')!, 130, 'Uniform');
      }, (e) => e.message.includes('not above the training rows (120)'));
    }]]},
  {name: 'Naive Bayes', fit: FITS.naiveBayes, ...NON_NATIVE, ...CLASSIFIER_ONLY, applicable: isNaiveBayesApplicable,
    interactive: isNaiveBayesInteractive, params: {classification: {model: 'gaussian_nb'}},
    shapes: {applicable: CLASSIFIER_APPLICABLE,
      interactive: [[stringShape(100000, 1), true], [stringShape(100001, 1), false]]}},
  {name: 'LDA', fit: FITS.lda, ...NON_NATIVE, ...CLASSIFIER_ONLY, applicable: isLdaApplicable,
    interactive: isLdaInteractive,
    params: {classification: {model: 'lda'}},
    shapes: {applicable: CLASSIFIER_APPLICABLE,
      interactive: [[stringShape(100000, 100), true], [stringShape(100001, 1), false],
        [stringShape(1000, 101), false]]},
    tests: [['Shrinkage', async () => {
      const {df, species} = await irisSpecies();
      expect((await paramsOf(trainLda(df, species, 0.1))).shrinkage, 0.1);
      expect('shrinkage' in await paramsOf(trainLda(df, species, 0)), false, 'shrinkage 0 must not be sent');
    }]]},
  {name: 'Neural Network', fit: FITS.mlp, ...NON_NATIVE, standard: true, seeded: true,
    floors: {species: 0.9, petal: 0.7}, applicable: isMlpApplicable, interactive: isMlpInteractive,
    params: both({hidden_layer_sizes: '100', activation: 'relu', max_iter: 200, learning_rate: 0.001, alpha: 0.0001,
      tol: 0.0001, random_state: 0, model: 'mlp'}),
    shapes: {applicable: BOOSTING_APPLICABLE,
      interactive: [[floatShape(2000, 50), true], [floatShape(2001, 1), false], [floatShape(1000, 51), false]]},
    tests: [['Hidden layers text', async () => {
      const {df, species} = await irisSpecies();
      const mlp = (layers: string) => trainMlp(df, species, layers, 'ReLU', 50, 0.001, 0.0001, 0);
      await expectRejected([['abc', () => mlp('abc')], ['0', () => mlp('0')], ['64,,32', () => mlp('64,,32')]]);
      expect((await paramsOf(mlp(' 64, 32 '))).hidden_layer_sizes, '64,32');
    }]]},
  {name: 'RealMLP', fit: FITS.realMlp, ...NON_NATIVE, seeded: true, floors: {species: 0.9, petal: 0.9},
    boolRows: 1000, applicable: isRealMlpApplicable, interactive: isRealMlpInteractive,
    params: both({model: 'realmlp', variant: 'td', n_epochs: 128, hidden_sizes: [256, 256, 256], batch_size: 256,
      val_fraction: 0.2, random_state: 0}),
    shapes: {applicable: [[floatShape(10000, 10), true], [floatShape(10001, 1), false], [floatShape(5000, 21), false]],
      interactive: [[floatShape(10, 1), false], [floatShape(150, 4), false]]},
    tests: [['Hidden layers text', async () => {
      const {df, species} = await irisSpecies();
      const realMlp = (layers: string) => trainRealMlp(df, species, 'TD', 4, layers, 0, 256, 0.2, 0);
      await expectRejected([['abc', () => realMlp('abc')], ['0', () => realMlp('0')],
        ['64,,32', () => realMlp('64,,32')]]);
      expect(JSON.stringify((await paramsOf(realMlp(' 128, 64 '))).hidden_sizes), JSON.stringify([128, 64]));
    }], ['Validation', async () => {
      const {df, species} = await irisSpecies();
      await expectRejected([
        ['epochs 0', () => trainRealMlp(df, species, 'TD', 0, '256,256,256', 0, 256, 0.2, 0)],
        ['batch size 0', () => trainRealMlp(df, species, 'TD', 4, '256,256,256', 0, 0, 0.2, 0)],
        ['learning rate -1', () => trainRealMlp(df, species, 'TD', 4, '256,256,256', -1, 256, 0.2, 0)],
        ['validation fraction 0', () => trainRealMlp(df, species, 'TD', 4, '256,256,256', 0, 256, 0, 0)],
        ['validation fraction 0.6', () => trainRealMlp(df, species, 'TD', 4, '256,256,256', 0, 256, 0.6, 0)],
      ]);
    }], ['Learning rate and variant', async () => {
      const {df, species} = await irisSpecies();
      const realMlp = (variant: string, lr: number) => trainRealMlp(df, species, variant, 4, '256,256,256', lr, 256,
        0.2, 0);
      expect('lr' in await paramsOf(realMlp('TD', 0)), false, 'learning rate 0 must not be sent');
      expect((await paramsOf(realMlp('TD', 0.05))).lr, 0.05);
      expect((await paramsOf(realMlp('TD-S', 0))).variant, 'td_s');
    }]]},
  {name: 'TabM', fit: FITS.tabm, ...TREES, floors: {species: 0.85, petal: 0.9},
    applicable: isTabmApplicable, interactive: isTabmInteractive,
    params: both({model: 'tabm', k: 8, d_block: 128, embeddings: 'piecewise-linear', max_epochs: 32, patience: 16,
      lr: 0.002, dropout: 0.1, val_fraction: 0.2, random_state: 0}),
    shapes: {applicable: [[floatShape(20000, 10), true], [floatShape(20001, 1), false], [floatShape(10000, 21), false]],
      interactive: [[floatShape(200, 10), true], [floatShape(201, 1), false], [floatShape(100, 21), false]]},
    tests: [['Validation', async () => {
      const {df, species} = await irisSpecies();
      const tabm = (k: number, width: number, epochs: number, patience: number, lr: number, dropout: number,
        fraction: number) => trainTabm(df, species, k, width, 'Piecewise linear', epochs, patience, lr, dropout,
        fraction, 0);
      await expectRejected([
        ['ensemble size 0', () => tabm(0, 128, 2, 16, 0.002, 0.1, 0.2)],
        ['width 0', () => tabm(8, 0, 2, 16, 0.002, 0.1, 0.2)],
        ['epochs 0', () => tabm(8, 128, 0, 16, 0.002, 0.1, 0.2)],
        ['patience 0', () => tabm(8, 128, 2, 0, 0.002, 0.1, 0.2)],
        ['learning rate 0', () => tabm(8, 128, 2, 16, 0, 0.1, 0.2)],
        ['dropout 1', () => tabm(8, 128, 2, 16, 0.002, 1, 0.2)],
        ['validation fraction 0', () => tabm(8, 128, 2, 16, 0.002, 0.1, 0)],
      ]);
    }], ['Embeddings', async () => {
      const {df, species} = await irisSpecies();
      const tabm = (embeddings: string) => trainTabm(df, species, 8, 128, embeddings, 2, 16, 0.002, 0.1, 0.2, 0);
      expect((await paramsOf(tabm('None'))).embeddings, 'none');
      expect((await paramsOf(tabm('Periodic'))).embeddings, 'periodic');
      expect((await paramsOf(tabm('Linear-ReLU'))).embeddings, 'linear-relu');
    }]]},
];

const sortedJson = (o: object) => JSON.stringify(Object.entries(o).sort(([a], [b]) => a.localeCompare(b)));

/** `features` zero FLOAT columns and a FLOAT / STRING target of `rows` rows: the predicates read only the shape */
function shapeOf(shape: Shape): {features: DG.ColumnList, target: DG.Column} {
  const cols = Array.from({length: shape.features}, (_, j) => DG.Column.fromType(DG.COLUMN_TYPE.FLOAT, `f${j}`,
    shape.rows));
  const type = shape.target === 'float' ? DG.COLUMN_TYPE.FLOAT : DG.COLUMN_TYPE.STRING;
  return {features: asTable(cols), target: DG.Column.fromType(type, 'target', shape.rows)};
}

function expectShapes(predicate: Predicate, shapes: [Shape, boolean][], label: string): void {
  for (const [shape, expected] of shapes) {
    const {features, target} = shapeOf(shape);
    expect(predicate(features, target), expected,
      `${label} at ${shape.rows} rows x ${shape.features} features (${shape.target} target)`);
  }
}

for (const m of METHODS) {
  category(m.name, () => {
    if (m.classification) {
      test('Iris classification round trip', async () => {
        const iris = await readIris();
        await roundTrip(m.name, m.fit, iris.col('Species')!, numericFeatures(iris, 'Species', ['col 1']),
          m.floors.species);
      }, {timeout: TIMEOUT});
    }

    if (m.regression) {
      test('Iris regression round trip', async () => {
        const iris = await readIris();
        await roundTrip(m.name, m.fit, iris.col('Petal.Length')!, numericFeatures(iris, 'Petal.Length', ['col 1']),
          m.floors.petal);
      }, {timeout: TIMEOUT});
    }

    test('Hyperparameters in the manifest', async () => {
      const {iris, df, species} = await irisSpecies();
      const targets: [string, DG.Column, object | undefined][] = [['classification', species, m.params.classification],
        ['regression', iris.col('Petal.Length')!, m.params.regression]];
      for (const [task, target, expected] of targets) {
        if (expected === undefined)
          continue;
        const manifest = readCruxManifest(await m.fit(df, target));
        expect(sortedJson(manifest.params), sortedJson(expected), `${task} params`);
        expect(JSON.stringify((manifest.pipeline ?? []).map((p) => p.step)),
          JSON.stringify(m.standard ? ['standard_scaler'] : []), `${task} pipeline`);
        expect(manifest.n_features, 4, `${task} features`);
      }
    }, {timeout: TIMEOUT});

    test('Missing feature values', async () => {
      const {iris, df, species} = await irisSpecies();
      if (m.nativeNans) {
        await expectGapsKept(m.fit, iris);
        return;
      }
      const gapped = irisMeasurements(iris, GAP_ROWS);
      await expectExceptionAsync(async () => {
        await m.fit(gapped, species);
      }, (e) => e.message.includes('missing feature values'));
      const prediction = await predictColumn(gapped, await m.fit(df, species));
      expect(prediction.stats.missingValueCount, GAP_ROWS.length, 'only the rows with gaps stay empty');
      for (const row of GAP_ROWS)
        expect(prediction.isNone(row), true, `row ${row} has a gap and must stay empty`);
    }, {timeout: TIMEOUT});

    test('Applicability', async () => {
      const iris = await readIris();
      const features = numericFeatures(iris, 'Species', ['col 1']);
      const species = iris.col('Species')!;
      const petal = iris.col('Petal.Length')!;
      expect(m.applicable(features, species), m.classification, 'iris Species');
      expect(m.applicable(numericFeatures(iris, 'Petal.Length', ['col 1']), petal), m.regression, 'iris Petal.Length');

      const target = m.classification ? species : petal;
      const withString = asTable([...features.toList().map(copyCol), copyCol(species)]);
      expect(m.applicable(withString, target), false, 'a STRING feature');
      const bigInt = DG.Column.fromBigInt64Array('big', new BigInt64Array(iris.rowCount));
      expect(m.applicable(features, bigInt), false, 'a BIG_INT target');
      if (!m.classification) {
        const bool = DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'bool', new Array(iris.rowCount).fill(true));
        expect(m.applicable(features, bool), false, 'a BOOL target');
      }

      expectShapes(m.applicable, m.shapes.applicable, 'applicable');
      expectShapes(m.interactive, m.shapes.interactive, 'interactive');
    }, {timeout: TIMEOUT});

    if (m.classification) {
      test('Bool target', async () => {
        const {df, target} = await demogFrame('CONTROL', ['AGE', 'HEIGHT', 'WEIGHT'], true);
        const mask = DG.BitSet.create(df.rowCount, (i) => i < (m.boolRows ?? df.rowCount));
        const frame = df.clone(mask);
        const bool = target.clone(mask);
        expect(bool.type, DG.COLUMN_TYPE.BOOL, 'CONTROL is expected to be a boolean column');
        const prediction = await applyTwice(frame, await m.fit(frame, bool), m.name);
        expect(prediction.type, DG.COLUMN_TYPE.BOOL, 'prediction of a bool target must be a BOOL column');
        expect(accuracy(bool, prediction) > 0.5, true, 'training accuracy is no better than chance');
      }, {timeout: TIMEOUT});
    }

    if (m.seeded) {
      test('Seed', async () => {
        const {iris, df, species} = await irisSpecies();
        const tasks: [DG.Column, boolean, number][] = [[species, m.classification, 0],
          [iris.col('Petal.Length')!, m.regression, 1e-6]];
        for (const [target, supported, tolerance] of tasks) {
          if (!supported)
            continue;
          const first = await predictColumn(df, await m.fit(df, target));
          const second = await predictColumn(df, await m.fit(df, target));
          expect(sameValues(first, second, tolerance), true, `${target.name}: the same seed gave another model`);
        }
      }, {timeout: TIMEOUT});
    }

    for (const [name, body] of m.tests ?? [])
      test(name, body, {timeout: TIMEOUT});
  });
}

type BenchCase = {name: string, features: DG.DataFrame, target: DG.Column, gaps: boolean};

let benchCases: Promise<BenchCase[]> | null = null;

/** A case of `CASES` in `baseline.ts`: its features and target */
async function baselineCase(id: string): Promise<BenchCase> {
  const {c, df, features} = await loadCase(id);
  return {name: id, features: df.clone(null, features), target: df.col(c.target)!, gaps: false};
}

function synthCase(name: string, df: DG.DataFrame, target: string): BenchCase {
  return {name, features: df, target: splitTarget(df, target), gaps: false};
}

async function loadBenchCases(): Promise<BenchCase[]> {
  const files = await Promise.all(['cars-price', 'iris-species', 'iris-petal', 'wine-type', 'wine-quality']
    .map(baselineCase));
  const demog = await demogFrame('SEX', ['AGE', 'HEIGHT', 'WEIGHT'], false);
  const cases = [...files, {name: 'demog-gaps', features: demog.df, target: demog.target, gaps: true}];
  if (!DG.Test.isInBenchmark)
    return cases;
  return [...cases, ...await Promise.all(['synth-cls-100k', 'synth-reg-100k'].map(baselineCase)),
    synthCase('synth-cls-20k-100', classificationDataset(20000, 100, true, 42), 'Labels'),
    synthCase('synth-reg-20k-100', seededRegression(20000, 100, 42), 'y'),
  ];
}

/** Absorbs the worker start (the first call, or the restart after the live-model cap) before a timed fit */
const warmUp = (iris: BenchCase) => FITS.naiveBayes(iris.features, iris.target);

/** Model size always; fit time of the cases the method calls interactive only in a benchmark run */
async function benchmark(m: MethodCase): Promise<void> {
  benchCases ??= loadBenchCases();
  const failures: string[] = [];
  const cases = await benchCases;
  const iris = cases.find((c) => c.name === 'iris-species')!;
  for (const c of cases) {
    if (c.gaps && !m.nativeNans)
      continue;
    const features = c.features.columns;
    if (!m.applicable(features, c.target)) {
      console.log(`BENCH ${m.name} ${c.name} skipped (not applicable)`);
      continue;
    }
    await warmUp(iris);
    let t0 = performance.now();
    const bytes = await m.fit(c.features, c.target);
    const fitMs = performance.now() - t0;
    t0 = performance.now();
    await predictColumn(c.features, bytes);
    const predictMs = performance.now() - t0;
    const interactive = m.interactive(features, c.target);
    console.log(`BENCH ${m.name} ${c.name} rows=${c.target.length} features=${features.length} ` +
      `fit=${fitMs.toFixed(0)} predict=${predictMs.toFixed(0)} blob=${bytes.length} interactive=${interactive}`);
    if (bytes.length > MAX_BLOB_BYTES)
      failures.push(`${c.name}: model of ${bytes.length} B above ${MAX_BLOB_BYTES} B`);
    if (DG.Test.isInBenchmark && interactive && fitMs > MAX_INTERACTIVE_FIT_MS)
      failures.push(`${c.name}: interactive, but the fit took ${fitMs.toFixed(0)} ms`);
  }
  expect(failures.length, 0, `${m.name}: ${failures.join('; ')}`);
}

category('Crux methods benchmark', () => {
  for (const m of METHODS)
    test(m.name, async () => benchmark(m), {timeout: BENCH_TIMEOUT, benchmark: true});
});
