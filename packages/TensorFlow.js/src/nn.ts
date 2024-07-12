import * as tf from '@tensorflow/tfjs';
import * as DG from 'datagrok-api/dg';
import { activations, processFeatures, processLables, saveModel } from './utils';
import { Initializer } from '@tensorflow/tfjs-layers/dist/initializers';
import { LayerNormalization, BatchNormalization} from '@tensorflow/tfjs-layers/dist/layers/normalization';
import { Regularizer } from '@tensorflow/tfjs-layers/dist/regularizers';

type RegularizationOptions = {
  kernelRegularizer?: Regularizer,
  biasRegularizer?: Regularizer,
  activityRegularizer?: Regularizer
}

export async function prepareAndTrainNN(
    dataframe: DG.DataFrame,
    predictColumn: DG.Column,
    activationFunc: activations = 'relu',
    hiddenLayersSizes: string = '128,128',
    initializer: string = 'randomUniform',
    constantValue: number = 0.5,
    normalMean: number = 0.5,
    normalStddev: number = 0.1,
    minValue: number = 0.0,
    maxValue: number = 1.0,
    varianceScale: number = 0.5,
    varianceFanningMode: "fanIn" | "fanOut" | "fanAvg" = "fanIn",
    varianceDistribution: 'normal'|'uniform'|'truncatedNormal' = "normal",
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
  // //TODO: Do it all with Danfo.js unless there is a better DataFrame interface
  // let trainX_df: DG.DataFrame = DG.DataFrame.fromColumns([]);
  // let trainY_df: DG.Column;
  // for (const [key, value] of Object.entries(data)) {
  //   if (key === predict_column) {
  //     trainY_df = DG.Column.fromList(value.type, key, value.data);
  //     continue
  //   }
  //   trainX_df.columns.add(DG.Column.fromList(value.type, key, value.data));
  // }

  // Encode categorical data
  const trainX_df = dataframe;
  const trainY_df = predictColumn;

  let trainX: tf.Tensor2D = processFeatures(trainX_df);
  let trainY: tf.Tensor2D = processLables(trainY_df!);

  const hiddenLayersSizesPrep = hiddenLayersSizes.split(",").map(val => Number.parseInt(val));
  
  //orthogonal and identity are not present
  let initializerPrep;
  switch (initializer) {
    case "constant":
      initializerPrep = tf.initializers.constant({value: constantValue});
      break;
    case "glorotNormal":
      initializerPrep = tf.initializers.glorotNormal({seed: seed});
      break;
    case "glorotUniform":
      initializerPrep = tf.initializers.glorotUniform({seed: seed});
      break;
    case "heNormal":
      initializerPrep = tf.initializers.heNormal({seed: seed});
      break;
    case "heUniform":
      initializerPrep = tf.initializers.heUniform({seed: seed});
      break;
    case "leCunNormal":
      initializerPrep = tf.initializers.leCunNormal({seed: seed});
      break;
    case "leCunUniform":
      initializerPrep = tf.initializers.leCunUniform({seed: seed});
      break;
    case "ones":
      initializerPrep = tf.initializers.ones();
      break;
    case "randomNormal":
      initializerPrep = tf.initializers.randomNormal({seed: seed, mean: normalMean, stddev: normalStddev});
      break;
    case "truncatedNormal":
      initializerPrep = tf.initializers.truncatedNormal({seed: seed, mean: normalMean, stddev: normalStddev});
      break;
    case "varianceScaling":
      initializerPrep = tf.initializers.varianceScaling({
        seed: seed,
        scale: varianceScale,
        mode: varianceFanningMode,
        distribution: varianceDistribution
      });
      break;
    case "zeros":
      initializerPrep = tf.initializers.zeros();
      break;
    //default is randomUniform
    case "randomUniform":
    default:
      initializerPrep = tf.initializers.randomUniform({seed: seed, minval: minValue, maxval: maxValue});
      break;
  }

  let innerNormPrep;
  switch (hiddenRepresentationNormalization) {
    case "layerNorm":
      innerNormPrep = tf.layers.layerNormalization({
        axis: normalizationAxis,
        epsilon: normalizationEpsilon,
        center: normalizationCenter,
        scale: normalizationScale
      });
      break;
    case "batchNorm":
      innerNormPrep = tf.layers.batchNormalization({
        axis: normalizationAxis,
        momentum: normalizationMomentum,
        epsilon: normalizationEpsilon,
        center: normalizationCenter,
        scale: normalizationScale
      });
      break;
    default:
      break;
  }

  let optimizerPrep;
  switch (optimizer) {
    case "Adadelta":
      optimizerPrep = tf.train.adadelta(learningRate=learningRate);
      break;
    case "Adagrad":
      optimizerPrep = tf.train.adagrad(learningRate=learningRate);
      break;
    case "Adamax":
      optimizerPrep = tf.train.adamax(learningRate=learningRate);
      break;
    //FIXME: need to specify momentum as well
    case "Momentum":
      optimizerPrep = tf.train.momentum(learningRate=learningRate);
      break;
    case "RMSProp":
      optimizerPrep = tf.train.rmsprop(learningRate=learningRate);
      break;
    case "SGD":
      optimizerPrep = tf.train.sgd(learningRate=learningRate);
      break;
    case "Adam":
    default:
      optimizerPrep = tf.train.adam(learningRate=learningRate);
      break;
  }

  let regularizationPrep: RegularizationOptions = {};
  if (kernelRegularization) {
    regularizationPrep.kernelRegularizer = tf.regularizers.l1l2({l1: regularizationL1, l2: regularizationL2});
  }
  if (biasRegularization) {
    regularizationPrep.biasRegularizer = tf.regularizers.l1l2({l1: regularizationL1, l2: regularizationL2});
  }
  if (activityRegularization) {
    regularizationPrep.activityRegularizer = tf.regularizers.l1l2({l1: regularizationL1, l2: regularizationL2});
  }

  const lossFuncPrep = lossFunc.split(",");
  const metricsPrep = metrics.split(",");

  const earlyStopping = {
    monitor: monitor,
    minDelta: minDelta,
    patience: patience,
    mode: mode,
    baseline: baseline
  };

  return await trainNeuralNetwork(
    trainX,
    trainY,
    activationFunc,
    hiddenLayersSizesPrep,
    initializerPrep,
    innerNormPrep,
    droputRate,
    regularizationPrep,
    finalActivation,
    optimizerPrep,
    lossFuncPrep,
    metricsPrep,
    batchSize,
    epochs,
    validationSplit,
    seed,
    earlyStopping
  );
}

export async function trainNeuralNetwork(
    X: tf.Tensor2D,
    Y: tf.Tensor2D,
    activationFunc: activations,
    hiddenLayersSizes: number[],
    initializer: Initializer,
    hiddenRepresentationNormalization: LayerNormalization | BatchNormalization | undefined,
    droputRate: number = 0.2,
    regularizationOptions: RegularizationOptions,
    finalActivation: activations,
    optimizer: tf.Optimizer,
    lossFunc: string | string[],
    metrics: string | string[],
    batchSize: number,
    epochs: number,
    validationSplit: number,
    seed: number,
    earlyStopping?: {
      monitor?: string,
      minDelta?: number,
      patience?: number,
      mode?: 'auto'|'min'|'max',
      baseline?: number,
      restoreBestWeights?: boolean
    }
  ) { 
  let model = tf.sequential();
  //FIXME: maybe there is a better way to obtain input shape?
  model.add(tf.layers.dense({inputShape: [X.shape[1]], units: hiddenLayersSizes[0], activation: activationFunc}));
  for (let layerUnits of hiddenLayersSizes.slice(1, hiddenLayersSizes.length)) {
      model.add(tf.layers.dense(Object.assign({
        units: layerUnits,
        activation: activationFunc,
        //TODO: seed must affect initializers
        kernelInitializer: initializer,
        biasInitializer: initializer
      }, regularizationOptions)));
      if (typeof hiddenRepresentationNormalization !== 'undefined') {
        model.add(hiddenRepresentationNormalization)
      }
      if (droputRate !== 0) {
        model.add(tf.layers.dropout({rate: droputRate, seed: seed}));
      }
  }
  model.add(tf.layers.dense({units: Y.shape[1]!, activation: finalActivation}));

  model.compile({optimizer: optimizer, loss: lossFunc, metrics: metrics});

  await model.fit(
    X,
    Y,
    {
      batchSize: batchSize,
      epochs: epochs,
      validationSplit: validationSplit,
      callbacks: tf.callbacks.earlyStopping(earlyStopping)
    }
  );

  return await saveModel(model);
}