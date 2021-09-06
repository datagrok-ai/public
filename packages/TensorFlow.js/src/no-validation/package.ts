import * as DG from 'datagrok-api/dg';

import * as tf from '@tensorflow/tfjs';
import { processFeatures, activations, loadModel } from './utils';
import { prepareAndTrainNN } from './nn';

export let _package = new DG.Package();

//name: trainNN
//description: Train a neural network using TensorFlow.js
//meta.mlname: tfjsNN
//meta.mlrole: train
//input: dataframe df {validators: containsMissingValues}
//input: string predict_column
//input: string activationFunc=relu {category: Parameters; choices: ["relu", "linear", "elu", "relu6", "selu", "sigmoid", "softplus", "softmax", "tanh", "softplus", "swish", "mish", "hardSigmoid"]} [Activation function to use.]
//input: string hiddenLayerSizes=16,32 {category: Parameters} [Number of units in each hidden layer.]
//input: string initializer=randomUniform {category: Initializer; choices: ["constant", "glorotNormal", "glorotUniform", "heNormal", "heUniform", "leCunNormal", "leCunUniform", "ones", "orthogonal", "randomNormal", "randomUniform", "truncatedNormal", "varianceScaling", "zeros"]} [Kernel and bias weights initializer method.]
//input: double constantValue=0.5 {category: Initializer} [Constant: The value for each element in the variable.]
//input: double normalMean=0.5 {category: Initializer} [Random & Truancated Normal: Mean of the random values to generate.]
//input: double normalStddev=0.1 {category: Initializer} [Random & Truancated Normal: Standard deviation of the random values to generate.]
//input: double minValue=0.0 {category: Initializer} [Random Uniform: Lower bound of the range of random values to generate.]
//input: double maxValue=1.0 {category: Initializer} [Random Uniform: Upper bound of the range of random values to generate.]
//input: double varianceScale=0.5 {category: Initializaer} [Variance Scaling: Scaling factor (positive float).]
//input: string varianceFanningMode=fanIn {category: Initializer; choices: ["fanIn", "fanOut", "fanAvg"]} [Variance Scaling: Fanning mode for inputs and outputs.]
//input: string varianceDistribution=normal {category: Initializer; choices: ["nomral", "uniform", "truncatedNormal"]} [Variance Scaling: Probabilistic distribution of the values.]
//input: string normalization=batchNorm {category: LayersNormalization; choices: ["none", "batchNorm", "layerNorm"]} [Normalization method to use.]
//input: int axis=-1 {category: LayersNormalization} [The integer axis that should be normalized (typically the features axis).]
//input: double momentum=0.99 {category: LayersNormalization} [Momentum of the moving average.]
//input: double epsilon=0.001 {category: LayersNormalization} [Small float added to the variance to avoid dividing by zero.]
//input: bool center=true {category: LayersNormalization} [If true, add offset of beta to normalized tensor. If false, beta is ignored.]
//input: bool scale=true {category: LayersNormalization} [ If true, multiply by gamma. If false, gamma is not used. When the next layer is linear (also e.g. nn.relu), this can be disabled since the scaling will be done by the next layer.]
//input: double droputRate=0.2 {category: Regularization} [Applies dropout to the layers. Float between 0 and 1. Fraction of the input units to drop. Set 0 to disable dropout.]
//input: double l1=0.01 {category: Regularization} [Adds a term to the loss to penalize large weights: loss += sum(l1 * abs(x))]
//input: double l2=0.01 {category: Regularization} [Adds a term to the loss to penalize large weights: loss += sum(l2 * x^2)]
//input: bool kernelRegularization=true {category: Regularization} [Apply regularizer function to the dense kernel weights matrix.]
//input: bool biasRegularization=true {category: Regularization} [Apply regularizer function to the bias vector.]
//input: bool activityRegularization=false {category: Regularization} [Apply regularizer function to the activation.]
//input: string finalActivation=relu {category: Parameters; choices: ["relu", "linear", "elu", "relu6", "selu", "sigmoid", "softplus", "softmax", "tanh", "softplus", "swish", "mish", "hardSigmoid"]} [Activation function to use for the output layer.]
//input: string optimizer=Adam {category: Parameters; choices: ['Adadelta', 'Adagrad', 'Adam', 'Adamax', 'Momentum', 'RMSProp', 'SGD']}
//input: double learningRate=0.001 {category: Parameters}
//input: string lossFunc=meanSquaredError {category: Parameters; choices: ['meanSquaredError', 'meanAbsoluteError', 'meanAbsolutePercentageError', 'meanSquaredLogarithmicError', 'squaredHinge', 'hinge', 'categoricalHinge', 'logcosh', 'categoricalCrossentropy', 'sparseCategoricalCrossentropy', 'binaryCrossentropy', 'kullbackLeiblerDivergence', 'poisson', 'cosineProximity']}
//input: string metrics=MSE {category: Parameters; choices: ['binaryAccuracy', 'categoricalAccuracy', 'precision', 'categoricalCrossentropy', 'sparseCategoricalCrossentropy', 'MSE', 'MAE', 'MAPE', 'cosine']}
//input: int batchSize=1 {category: Parameters}
//input: int epochs=10 {category: Parameters}
//input: double validationSplit=0.2 {category: Parameters}
//input: int seed=42
//input: string constraint {choices: ["maxNorm", "minMaxNorm", "nonNeg", "unitNorm", ]}
//input: string monitor=val_loss {category: EarlyStopping}
//input: double minDelta=0.001 {category: EarlyStopping}
//input: int patience=5 {category: EarlyStopping}
//input: string mode=auto {category: EarlyStopping; choices: ["auto", "min", "max"]}
//input: double baseline=0.5 {category: EarlyStopping}
//output: dynamic model
export async function trainNN(
    df: DG.DataFrame,
    predict_column: string,
    activationFunc: activations = 'relu',
    hiddenLayersSizes: string = '128,128',
    initializer: string = 'randomUniform',
    constantValue: number = 0.5,
    normalMean: number = 0.5,
    normalStddev: number = 0.1,
    minValue: number = 0.0,
    maxValue: number = 1.0,
    varianceScale: number = 0.5,
    varianceFanningMode: string = "fanIn",
    varianceDistribution: string = "normal",
    hiddenRepresentationNormalization: string = "batchNorm",
    normalizationAxis: number = -1,
    normalizationMomentum: number = 0.99,
    normalizationEpsilon: number = 0.001,
    normalizationCenter: boolean = true,
    normalizationScale: boolean = true,
    droputRate: number = 0.2,
    regularizationL1: number = 0.01,
    regularizationL2: number = 0.01,
    kernelRegularization: boolean = true,
    biasRegularization: boolean = true,
    activityRegularization: boolean = false,
    finalActivation: activations = 'relu',
    optimizer: string = 'Adam',
    learningRate: number = 0.001,
    lossFunc: string = 'meanSquaredError',
    metrics: string = 'MSE',
    batchSize: number = 1,
    epochs: number = 1,
    validationSplit: number,
    seed: number = 42,
    monitor?: string,
    minDelta?: number,
    patience?: number,
    mode?: 'auto'|'min'|'max',
    baseline?: number
) {
  //TODO: reorganize
  // let columnsObject: object = {};
  // df.columns.toList().forEach((col: DG.Column) => {columnsObject[col.name] = {type: col.type, data: col.toList() } });
  let result = await prepareAndTrainNN(
    df,
    predict_column,
    activationFunc,
    hiddenLayersSizes,
    initializer,
    constantValue,
    normalMean,
    normalStddev,
    minValue,
    maxValue,
    varianceScale,
    varianceFanningMode,
    varianceDistribution,
    hiddenRepresentationNormalization,
    normalizationAxis,
    normalizationMomentum,
    normalizationEpsilon,
    normalizationCenter,
    normalizationScale,
    droputRate,
    regularizationL1,
    regularizationL2,
    kernelRegularization,
    biasRegularization,
    activityRegularization,
    finalActivation,
    optimizer,
    learningRate,
    lossFunc,
    metrics,
    batchSize,
    epochs,
    validationSplit,
    seed,
    monitor,
    minDelta,
    patience,
    mode,
    baseline
  );

  return result;
}

//name: applyNN
//meta.mlname: tfjsNN
//meta.mlrole: apply
//input: dataframe df
//input: dynamic model
//output: dataframe data_out
export async function applyNN(
  df: DG.DataFrame,
  model: Uint8Array
) {
  let loadedModel = await loadModel(model);
  let testX = processFeatures(df);
  let predictions: any[];
  // if (prediction is not category) {
  predictions = loadedModel.predict(testX).arraySync().flat();
  // } else { // is category
  // predictions = predictions.map((vector) => oneHotToCategory(vector, categoryArray));
  // }

  return DG.DataFrame.fromColumns([DG.Column.fromList(DG.TYPE.FLOAT, "pred", predictions)]);
}
