import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Scripts {

}

export namespace Funcs {
  export async function info(): Promise<any> {
    return await grok.functions.call('Eda:Info', {});
  }

  export async function init(): Promise<any> {
    return await grok.functions.call('Eda:Init', {});
  }

  //Density-based spatial clustering of applications with noise (DBSCAN)
  export async function dbScan(df: DG.DataFrame, xCol: DG.Column, yCol: DG.Column, epsilon: number, minPts: number): Promise<any> {
    return await grok.functions.call('Eda:DbScan', { df, xCol, yCol, epsilon, minPts });
  }

  //Principal component analysis (PCA)
  export async function pca(table: DG.DataFrame, features: string[], components: number, center: boolean, scale: boolean): Promise<any> {
    return await grok.functions.call('Eda:PCA', { table, features, components, center, scale });
  }

  export async function dbscanPostProcessingFunction(col1: DG.Column, col2: DG.Column, epsilon: number, minimumPoints: number): Promise<any> {
    return await grok.functions.call('Eda:DbscanPostProcessingFunction', { col1, col2, epsilon, minimumPoints });
  }

  export async function numberPreprocessingFunction(col: DG.Column, _metric: string): Promise<any> {
    return await grok.functions.call('Eda:NumberPreprocessingFunction', { col, _metric });
  }

  export async function stringPreprocessingFunction(col: DG.Column, _metric: string): Promise<any> {
    return await grok.functions.call('Eda:StringPreprocessingFunction', { col, _metric });
  }

  export async function reduceDimensionality(): Promise<any> {
    return await grok.functions.call('Eda:ReduceDimensionality', {});
  }

  export async function getMCLEditor(call: any): Promise<any> {
    return await grok.functions.call('Eda:GetMCLEditor', { call });
  }

  //Markov clustering (MCL) is an unsupervised clustering algorithm for graphs based on simulation of stochastic flow.
  export async function mclclustering(df: DG.DataFrame, aggregationMethod: string, preprocessingFuncArgs: any, threshold: number, maxIterations: number, useWebGPU: boolean, inflate: number, minClusterSize: number): Promise<any> {
    return await grok.functions.call('Eda:MCLClustering', { df, aggregationMethod, preprocessingFuncArgs, threshold, maxIterations, useWebGPU, inflate, minClusterSize });
  }

  //Markov clustering viewer
  export async function markovClusteringViewer(): Promise<any> {
    return await grok.functions.call('Eda:MarkovClusteringViewer', {});
  }

  //Compute partial least squares (PLS) regression analysis components: prediction, regression coefficients, T- & U-scores, X-loadings.
  export async function pls(table: DG.DataFrame, features: string[], predict: DG.Column, components: number, names: DG.Column): Promise<any> {
    return await grok.functions.call('Eda:PLS', { table, features, predict, components, names });
  }

  //Compute partial least squares (PLS) regression components. They maximally summarize the variation of the predictors while maximizing correlation with the response variable.
  export async function topMenuPLS(): Promise<any> {
    return await grok.functions.call('Eda:TopMenuPLS', {});
  }

  //Multidimensional data analysis using partial least squares (PLS) regression.
  export async function mva(): Promise<any> {
    return await grok.functions.call('Eda:MVA', {});
  }

  //Multidimensional data analysis using partial least squares (PLS) regression. It identifies latent factors and constructs a linear model based on them.
  export async function demoMultivariateAnalysis(): Promise<any> {
    return await grok.functions.call('Eda:DemoMultivariateAnalysis', {});
  }

  export async function trainLinearKernelSVM(df: DG.DataFrame, predictColumn: DG.Column, gamma: number): Promise<any> {
    return await grok.functions.call('Eda:TrainLinearKernelSVM', { df, predictColumn, gamma });
  }

  export async function applyLinearKernelSVM(df: DG.DataFrame, model: any): Promise<any> {
    return await grok.functions.call('Eda:ApplyLinearKernelSVM', { df, model });
  }

  export async function isApplicableLinearKernelSVM(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('Eda:IsApplicableLinearKernelSVM', { df, predictColumn });
  }

  export async function isInteractiveLinearKernelSVM(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('Eda:IsInteractiveLinearKernelSVM', { df, predictColumn });
  }

  export async function visualizeLinearKernelSVM(df: DG.DataFrame, targetColumn: DG.Column, predictColumn: DG.Column, model: any): Promise<any> {
    return await grok.functions.call('Eda:VisualizeLinearKernelSVM', { df, targetColumn, predictColumn, model });
  }

  export async function trainRBFkernelSVM(df: DG.DataFrame, predictColumn: DG.Column, gamma: number, sigma: number): Promise<any> {
    return await grok.functions.call('Eda:TrainRBFkernelSVM', { df, predictColumn, gamma, sigma });
  }

  export async function applyRBFkernelSVM(df: DG.DataFrame, model: any): Promise<any> {
    return await grok.functions.call('Eda:ApplyRBFkernelSVM', { df, model });
  }

  export async function isApplicableRBFkernelSVM(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('Eda:IsApplicableRBFkernelSVM', { df, predictColumn });
  }

  export async function isInteractiveRBFkernelSVM(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('Eda:IsInteractiveRBFkernelSVM', { df, predictColumn });
  }

  export async function visualizeRBFkernelSVM(df: DG.DataFrame, targetColumn: DG.Column, predictColumn: DG.Column, model: any): Promise<any> {
    return await grok.functions.call('Eda:VisualizeRBFkernelSVM', { df, targetColumn, predictColumn, model });
  }

  export async function trainPolynomialKernelSVM(df: DG.DataFrame, predictColumn: DG.Column, gamma: number, c: number, d: number): Promise<any> {
    return await grok.functions.call('Eda:TrainPolynomialKernelSVM', { df, predictColumn, gamma, c, d });
  }

  export async function applyPolynomialKernelSVM(df: DG.DataFrame, model: any): Promise<any> {
    return await grok.functions.call('Eda:ApplyPolynomialKernelSVM', { df, model });
  }

  export async function isApplicablePolynomialKernelSVM(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('Eda:IsApplicablePolynomialKernelSVM', { df, predictColumn });
  }

  export async function isInteractivePolynomialKernelSVM(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('Eda:IsInteractivePolynomialKernelSVM', { df, predictColumn });
  }

  export async function visualizePolynomialKernelSVM(df: DG.DataFrame, targetColumn: DG.Column, predictColumn: DG.Column, model: any): Promise<any> {
    return await grok.functions.call('Eda:VisualizePolynomialKernelSVM', { df, targetColumn, predictColumn, model });
  }

  export async function trainSigmoidKernelSVM(df: DG.DataFrame, predictColumn: DG.Column, gamma: number, kappa: number, theta: number): Promise<any> {
    return await grok.functions.call('Eda:TrainSigmoidKernelSVM', { df, predictColumn, gamma, kappa, theta });
  }

  export async function applySigmoidKernelSVM(df: DG.DataFrame, model: any): Promise<any> {
    return await grok.functions.call('Eda:ApplySigmoidKernelSVM', { df, model });
  }

  export async function isApplicableSigmoidKernelSVM(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('Eda:IsApplicableSigmoidKernelSVM', { df, predictColumn });
  }

  export async function isInteractiveSigmoidKernelSVM(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('Eda:IsInteractiveSigmoidKernelSVM', { df, predictColumn });
  }

  export async function visualizeSigmoidKernelSVM(df: DG.DataFrame, targetColumn: DG.Column, predictColumn: DG.Column, model: any): Promise<any> {
    return await grok.functions.call('Eda:VisualizeSigmoidKernelSVM', { df, targetColumn, predictColumn, model });
  }

  //One-way analysis of variances (ANOVA) determines whether the examined factor has a significant impact on the explored feature.
  export async function anova(): Promise<any> {
    return await grok.functions.call('Eda:Anova', {});
  }

  //Missing values imputation using the k-nearest neighbors method (KNN)
  export async function kNNImputation(): Promise<any> {
    return await grok.functions.call('Eda:KNNImputation', {});
  }

  //Missing values imputation using the k-nearest neighbors method
  export async function kNNImputationForTable(table: DG.DataFrame): Promise<any> {
    return await grok.functions.call('Eda:KNNImputationForTable', { table });
  }

  export async function trainLinearRegression(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('Eda:TrainLinearRegression', { df, predictColumn });
  }

  export async function applyLinearRegression(df: DG.DataFrame, model: any): Promise<any> {
    return await grok.functions.call('Eda:ApplyLinearRegression', { df, model });
  }

  export async function isApplicableLinearRegression(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('Eda:IsApplicableLinearRegression', { df, predictColumn });
  }

  export async function isInteractiveLinearRegression(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('Eda:IsInteractiveLinearRegression', { df, predictColumn });
  }

  export async function trainSoftmax(df: DG.DataFrame, predictColumn: DG.Column, rate: number, iterations: number, penalty: number, tolerance: number): Promise<any> {
    return await grok.functions.call('Eda:TrainSoftmax', { df, predictColumn, rate, iterations, penalty, tolerance });
  }

  export async function applySoftmax(df: DG.DataFrame, model: any): Promise<any> {
    return await grok.functions.call('Eda:ApplySoftmax', { df, model });
  }

  export async function isApplicableSoftmax(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('Eda:IsApplicableSoftmax', { df, predictColumn });
  }

  export async function isInteractiveSoftmax(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('Eda:IsInteractiveSoftmax', { df, predictColumn });
  }

  export async function trainPLSRegression(df: DG.DataFrame, predictColumn: DG.Column, components: number): Promise<any> {
    return await grok.functions.call('Eda:TrainPLSRegression', { df, predictColumn, components });
  }

  export async function applyPLSRegression(df: DG.DataFrame, model: any): Promise<any> {
    return await grok.functions.call('Eda:ApplyPLSRegression', { df, model });
  }

  export async function isApplicablePLSRegression(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('Eda:IsApplicablePLSRegression', { df, predictColumn });
  }

  export async function visualizePLSRegression(df: DG.DataFrame, targetColumn: DG.Column, predictColumn: DG.Column, model: any): Promise<any> {
    return await grok.functions.call('Eda:VisualizePLSRegression', { df, targetColumn, predictColumn, model });
  }

  export async function isInteractivePLSRegression(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('Eda:IsInteractivePLSRegression', { df, predictColumn });
  }

  export async function trainXGBooster(df: DG.DataFrame, predictColumn: DG.Column, iterations: number, eta: number, maxDepth: number, lambda: number, alpha: number): Promise<any> {
    return await grok.functions.call('Eda:TrainXGBooster', { df, predictColumn, iterations, eta, maxDepth, lambda, alpha });
  }

  export async function applyXGBooster(df: DG.DataFrame, model: any): Promise<any> {
    return await grok.functions.call('Eda:ApplyXGBooster', { df, model });
  }

  export async function isInteractiveXGBooster(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('Eda:IsInteractiveXGBooster', { df, predictColumn });
  }

  export async function isApplicableXGBooster(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('Eda:IsApplicableXGBooster', { df, predictColumn });
  }
}
