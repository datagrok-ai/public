import * as tf from '@tensorflow/tfjs';
import { processFeatures, processLables, saveModel } from './utils';
export async function prepareAndTrainNN(dataframe, predictColumn, activationFunc = 'relu', hiddenLayersSizes = '128,128', initializer = 'randomUniform', constantValue = 0.5, normalMean = 0.5, normalStddev = 0.1, minValue = 0.0, maxValue = 1.0, varianceScale = 0.5, varianceFanningMode = "fanIn", varianceDistribution = "normal", hiddenRepresentationNormalization = "batchNorm", normalizationAxis = -1, normalizationMomentum = 0.99, normalizationEpsilon = 0.001, normalizationCenter = true, normalizationScale = true, droputRate = 0.2, regularizationL1 = 0.01, regularizationL2 = 0.01, kernelRegularization = true, biasRegularization = true, activityRegularization = false, finalActivation = 'relu', optimizer = 'Adam', learningRate = 0.001, lossFunc = 'meanSquaredError', metrics = 'MSE', batchSize = 1, epochs = 1, validationSplit, seed = 42, monitor, minDelta, patience, mode, baseline) {
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
    let trainX = processFeatures(trainX_df);
    let trainY = processLables(trainY_df);
    const hiddenLayersSizesPrep = hiddenLayersSizes.split(",").map(val => Number.parseInt(val));
    //orthogonal and identity are not present
    let initializerPrep;
    switch (initializer) {
        case "constant":
            initializerPrep = tf.initializers.constant({ value: constantValue });
            break;
        case "glorotNormal":
            initializerPrep = tf.initializers.glorotNormal({ seed: seed });
            break;
        case "glorotUniform":
            initializerPrep = tf.initializers.glorotUniform({ seed: seed });
            break;
        case "heNormal":
            initializerPrep = tf.initializers.heNormal({ seed: seed });
            break;
        case "heUniform":
            initializerPrep = tf.initializers.heUniform({ seed: seed });
            break;
        case "leCunNormal":
            initializerPrep = tf.initializers.leCunNormal({ seed: seed });
            break;
        case "leCunUniform":
            initializerPrep = tf.initializers.leCunUniform({ seed: seed });
            break;
        case "ones":
            initializerPrep = tf.initializers.ones();
            break;
        case "randomNormal":
            initializerPrep = tf.initializers.randomNormal({ seed: seed, mean: normalMean, stddev: normalStddev });
            break;
        case "truncatedNormal":
            initializerPrep = tf.initializers.truncatedNormal({ seed: seed, mean: normalMean, stddev: normalStddev });
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
            initializerPrep = tf.initializers.randomUniform({ seed: seed, minval: minValue, maxval: maxValue });
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
        default:
            innerNormPrep = tf.layers.batchNormalization({
                axis: normalizationAxis,
                momentum: normalizationMomentum,
                epsilon: normalizationEpsilon,
                center: normalizationCenter,
                scale: normalizationScale
            });
            break;
    }
    let optimizerPrep;
    switch (optimizer) {
        case "Adadelta":
            optimizerPrep = tf.train.adadelta(learningRate = learningRate);
            break;
        case "Adagrad":
            optimizerPrep = tf.train.adagrad(learningRate = learningRate);
            break;
        case "Adamax":
            optimizerPrep = tf.train.adamax(learningRate = learningRate);
            break;
        //FIXME: need to specify momentum as well
        case "Momentum":
            optimizerPrep = tf.train.momentum(learningRate = learningRate);
            break;
        case "RMSProp":
            optimizerPrep = tf.train.rmsprop(learningRate = learningRate);
            break;
        case "SGD":
            optimizerPrep = tf.train.sgd(learningRate = learningRate);
            break;
        case "Adam":
        default:
            optimizerPrep = tf.train.adam(learningRate = learningRate);
            break;
    }
    let regularizationPrep = {};
    if (kernelRegularization) {
        regularizationPrep.kernelRegularizer = tf.regularizers.l1l2({ l1: regularizationL1, l2: regularizationL2 });
    }
    if (biasRegularization) {
        regularizationPrep.biasRegularizer = tf.regularizers.l1l2({ l1: regularizationL1, l2: regularizationL2 });
    }
    if (activityRegularization) {
        regularizationPrep.activityRegularizer = tf.regularizers.l1l2({ l1: regularizationL1, l2: regularizationL2 });
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
    return await trainNeuralNetwork(trainX, trainY, activationFunc, hiddenLayersSizesPrep, initializerPrep, innerNormPrep, droputRate, regularizationPrep, finalActivation, optimizerPrep, lossFuncPrep, metricsPrep, batchSize, epochs, validationSplit, seed, earlyStopping);
}
export async function trainNeuralNetwork(X, Y, activationFunc, hiddenLayersSizes, initializer, hiddenRepresentationNormalization, droputRate = 0.2, regularizationOptions, finalActivation, optimizer, lossFunc, metrics, batchSize, epochs, validationSplit, seed, earlyStopping) {
    let model = tf.sequential();
    //FIXME: maybe there is a better way to obtain input shape?
    model.add(tf.layers.dense({ inputShape: [X.shape[1]], units: hiddenLayersSizes[0], activation: activationFunc }));
    for (let layerUnits of hiddenLayersSizes.slice(1, hiddenLayersSizes.length)) {
        model.add(tf.layers.dense(Object.assign({
            units: layerUnits,
            activation: activationFunc,
            //TODO: seed must affect initializers
            kernelInitializer: initializer,
            biasInitializer: initializer
        }, regularizationOptions)));
        if (typeof hiddenRepresentationNormalization !== 'undefined') {
            model.add(hiddenRepresentationNormalization);
        }
        if (droputRate !== 0) {
            model.add(tf.layers.dropout({ rate: droputRate, seed: seed }));
        }
    }
    model.add(tf.layers.dense({ units: Y.shape[1], activation: finalActivation }));
    model.compile({ optimizer: optimizer, loss: lossFunc, metrics: metrics });
    await model.fit(X, Y, {
        batchSize: batchSize,
        epochs: epochs,
        validationSplit: validationSplit,
        callbacks: tf.callbacks.earlyStopping(earlyStopping)
    });
    return await saveModel(model);
}
