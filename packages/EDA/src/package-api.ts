import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace scripts {

}

export namespace funcs {
  export async function info(): Promise<any> {
    return await grok.functions.call('@datagrok/eda:Info', {});
  }

  export async function init(): Promise<any> {
    return await grok.functions.call('@datagrok/eda:Init', {});
  }

  //Density-based spatial clustering of applications with noise (DBSCAN)
  export async function dbScan(df: DG.DataFrame, xCol: DG.Column, yCol: DG.Column, epsilon: number, minPts: number): Promise<any> {
    return await grok.functions.call('@datagrok/eda:DbScan', { df, xCol, yCol, epsilon, minPts });
  }

  //Principal component analysis (PCA)
  export async function pca(table: DG.DataFrame, features: string[], components: number, center: boolean, scale: boolean): Promise<any> {
    return await grok.functions.call('@datagrok/eda:PCA', { table, features, components, center, scale });
  }

  export async function dbscanPostProcessingFunction(col1: DG.Column, col2: DG.Column, epsilon: number, minimumPoints: number): Promise<any> {
    return await grok.functions.call('@datagrok/eda:DbscanPostProcessingFunction', { col1, col2, epsilon, minimumPoints });
  }

  export async function numberPreprocessingFunction(col: DG.Column, _metric: string): Promise<any> {
    return await grok.functions.call('@datagrok/eda:NumberPreprocessingFunction', { col, _metric });
  }

  export async function stringPreprocessingFunction(col: DG.Column, _metric: string): Promise<any> {
    return await grok.functions.call('@datagrok/eda:StringPreprocessingFunction', { col, _metric });
  }

  export async function reduceDimensionality(): Promise<any> {
    return await grok.functions.call('@datagrok/eda:ReduceDimensionality', {});
  }

  export async function getMCLEditor(call: any): Promise<any> {
    return await grok.functions.call('@datagrok/eda:GetMCLEditor', { call });
  }

  //Markov clustering (MCL) is an unsupervised clustering algorithm for graphs based on simulation of stochastic flow.
  export async function mclclustering(df: DG.DataFrame, aggregationMethod: string, preprocessingFuncArgs: any, threshold: number, maxIterations: number, useWebGPU: boolean, inflate: number, minClusterSize: number): Promise<any> {
    return await grok.functions.call('@datagrok/eda:MCLClustering', { df, aggregationMethod, preprocessingFuncArgs, threshold, maxIterations, useWebGPU, inflate, minClusterSize });
  }

  //Markov clustering viewer
  export async function markovClusteringViewer(): Promise<any> {
    return await grok.functions.call('@datagrok/eda:MarkovClusteringViewer', {});
  }

  //Compute partial least squares (PLS) regression analysis components: prediction, regression coefficients, T- & U-scores, X-loadings.
  export async function pls(table: DG.DataFrame, features: string[], predict: DG.Column, components: number, names: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/eda:PLS', { table, features, predict, components, names });
  }

  //Compute partial least squares (PLS) regression components. They maximally summarize the variation of the predictors while maximizing correlation with the response variable.
  export async function topMenuPLS(): Promise<any> {
    return await grok.functions.call('@datagrok/eda:TopMenuPLS', {});
  }

  //Multidimensional data analysis using partial least squares (PLS) regression.
  export async function mva(): Promise<any> {
    return await grok.functions.call('@datagrok/eda:MVA', {});
  }

  //Multidimensional data analysis using partial least squares (PLS) regression. It identifies latent factors and constructs a linear model based on them.
  export async function demoMultivariateAnalysis(): Promise<any> {
    return await grok.functions.call('@datagrok/eda:DemoMultivariateAnalysis', {});
  }

  export async function trainLinearKernelSVM(df: DG.DataFrame, predictColumn: DG.Column, gamma: number): Promise<any> {
    return await grok.functions.call('@datagrok/eda:TrainLinearKernelSVM', { df, predictColumn, gamma });
  }

  export async function applyLinearKernelSVM(df: DG.DataFrame, model: any): Promise<any> {
    return await grok.functions.call('@datagrok/eda:ApplyLinearKernelSVM', { df, model });
  }

  export async function isApplicableLinearKernelSVM(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/eda:IsApplicableLinearKernelSVM', { df, predictColumn });
  }

  export async function isInteractiveLinearKernelSVM(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/eda:IsInteractiveLinearKernelSVM', { df, predictColumn });
  }

  export async function visualizeLinearKernelSVM(df: DG.DataFrame, targetColumn: DG.Column, predictColumn: DG.Column, model: any): Promise<any> {
    return await grok.functions.call('@datagrok/eda:VisualizeLinearKernelSVM', { df, targetColumn, predictColumn, model });
  }

  export async function trainRBFkernelSVM(df: DG.DataFrame, predictColumn: DG.Column, gamma: number, sigma: number): Promise<any> {
    return await grok.functions.call('@datagrok/eda:TrainRBFkernelSVM', { df, predictColumn, gamma, sigma });
  }

  export async function applyRBFkernelSVM(df: DG.DataFrame, model: any): Promise<any> {
    return await grok.functions.call('@datagrok/eda:ApplyRBFkernelSVM', { df, model });
  }

  export async function isApplicableRBFkernelSVM(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/eda:IsApplicableRBFkernelSVM', { df, predictColumn });
  }

  export async function isInteractiveRBFkernelSVM(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/eda:IsInteractiveRBFkernelSVM', { df, predictColumn });
  }

  export async function visualizeRBFkernelSVM(df: DG.DataFrame, targetColumn: DG.Column, predictColumn: DG.Column, model: any): Promise<any> {
    return await grok.functions.call('@datagrok/eda:VisualizeRBFkernelSVM', { df, targetColumn, predictColumn, model });
  }

  export async function trainPolynomialKernelSVM(df: DG.DataFrame, predictColumn: DG.Column, gamma: number, c: number, d: number): Promise<any> {
    return await grok.functions.call('@datagrok/eda:TrainPolynomialKernelSVM', { df, predictColumn, gamma, c, d });
  }

  export async function applyPolynomialKernelSVM(df: DG.DataFrame, model: any): Promise<any> {
    return await grok.functions.call('@datagrok/eda:ApplyPolynomialKernelSVM', { df, model });
  }

  export async function isApplicablePolynomialKernelSVM(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/eda:IsApplicablePolynomialKernelSVM', { df, predictColumn });
  }

  export async function isInteractivePolynomialKernelSVM(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/eda:IsInteractivePolynomialKernelSVM', { df, predictColumn });
  }

  export async function visualizePolynomialKernelSVM(df: DG.DataFrame, targetColumn: DG.Column, predictColumn: DG.Column, model: any): Promise<any> {
    return await grok.functions.call('@datagrok/eda:VisualizePolynomialKernelSVM', { df, targetColumn, predictColumn, model });
  }

  export async function trainSigmoidKernelSVM(df: DG.DataFrame, predictColumn: DG.Column, gamma: number, kappa: number, theta: number): Promise<any> {
    return await grok.functions.call('@datagrok/eda:TrainSigmoidKernelSVM', { df, predictColumn, gamma, kappa, theta });
  }

  export async function applySigmoidKernelSVM(df: DG.DataFrame, model: any): Promise<any> {
    return await grok.functions.call('@datagrok/eda:ApplySigmoidKernelSVM', { df, model });
  }

  export async function isApplicableSigmoidKernelSVM(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/eda:IsApplicableSigmoidKernelSVM', { df, predictColumn });
  }

  export async function isInteractiveSigmoidKernelSVM(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/eda:IsInteractiveSigmoidKernelSVM', { df, predictColumn });
  }

  export async function visualizeSigmoidKernelSVM(df: DG.DataFrame, targetColumn: DG.Column, predictColumn: DG.Column, model: any): Promise<any> {
    return await grok.functions.call('@datagrok/eda:VisualizeSigmoidKernelSVM', { df, targetColumn, predictColumn, model });
  }

  //One-way analysis of variances (ANOVA) determines whether the examined factor has a significant impact on the explored feature.
  export async function anova(): Promise<any> {
    return await grok.functions.call('@datagrok/eda:Anova', {});
  }

  //Missing values imputation using the k-nearest neighbors method (KNN)
  export async function kNNImputation(): Promise<any> {
    return await grok.functions.call('@datagrok/eda:KNNImputation', {});
  }

  //Missing values imputation using the k-nearest neighbors method
  export async function kNNImputationForTable(table: DG.DataFrame): Promise<any> {
    return await grok.functions.call('@datagrok/eda:KNNImputationForTable', { table });
  }

  export async function trainLinearRegression(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/eda:TrainLinearRegression', { df, predictColumn });
  }

  export async function applyLinearRegression(df: DG.DataFrame, model: any): Promise<any> {
    return await grok.functions.call('@datagrok/eda:ApplyLinearRegression', { df, model });
  }

  export async function isApplicableLinearRegression(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/eda:IsApplicableLinearRegression', { df, predictColumn });
  }

  export async function isInteractiveLinearRegression(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/eda:IsInteractiveLinearRegression', { df, predictColumn });
  }

  export async function trainSoftmax(df: DG.DataFrame, predictColumn: DG.Column, rate: number, iterations: number, penalty: number, tolerance: number): Promise<any> {
    return await grok.functions.call('@datagrok/eda:TrainSoftmax', { df, predictColumn, rate, iterations, penalty, tolerance });
  }

  export async function applySoftmax(df: DG.DataFrame, model: any): Promise<any> {
    return await grok.functions.call('@datagrok/eda:ApplySoftmax', { df, model });
  }

  export async function isApplicableSoftmax(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/eda:IsApplicableSoftmax', { df, predictColumn });
  }

  export async function isInteractiveSoftmax(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/eda:IsInteractiveSoftmax', { df, predictColumn });
  }

  export async function trainPLSRegression(df: DG.DataFrame, predictColumn: DG.Column, components: number): Promise<any> {
    return await grok.functions.call('@datagrok/eda:TrainPLSRegression', { df, predictColumn, components });
  }

  export async function applyPLSRegression(df: DG.DataFrame, model: any): Promise<any> {
    return await grok.functions.call('@datagrok/eda:ApplyPLSRegression', { df, model });
  }

  export async function isApplicablePLSRegression(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/eda:IsApplicablePLSRegression', { df, predictColumn });
  }

  export async function visualizePLSRegression(df: DG.DataFrame, targetColumn: DG.Column, predictColumn: DG.Column, model: any): Promise<any> {
    return await grok.functions.call('@datagrok/eda:VisualizePLSRegression', { df, targetColumn, predictColumn, model });
  }

  export async function isInteractivePLSRegression(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/eda:IsInteractivePLSRegression', { df, predictColumn });
  }

  export async function trainXGBooster(df: DG.DataFrame, predictColumn: DG.Column, iterations: number, eta: number, maxDepth: number, lambda: number, alpha: number): Promise<any> {
    return await grok.functions.call('@datagrok/eda:TrainXGBooster', { df, predictColumn, iterations, eta, maxDepth, lambda, alpha });
  }

  export async function applyXGBooster(df: DG.DataFrame, model: any): Promise<any> {
    return await grok.functions.call('@datagrok/eda:ApplyXGBooster', { df, model });
  }

  export async function isInteractiveXGBooster(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/eda:IsInteractiveXGBooster', { df, predictColumn });
  }

  export async function isApplicableXGBooster(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/eda:IsApplicableXGBooster', { df, predictColumn });
  }
}
