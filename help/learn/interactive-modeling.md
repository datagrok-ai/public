---
title: "Interactive modeling"
---

The predictive modeling toolkit allows an interactive visualization tool for the models. 

There are a lot of cases when models do not have a lot of parameters and could be trained relatively fast. For such situation, Datagrok builds a UI for the model: it has all the parameters that model uses. The user has the ability to change parameters, and the model will be retrained, providing the user with visualizations

![](./interactive-modeling.gif)

## Workflow

To start setting up a model, we start with the data configuration. First, we select a table to work with, then we have to choose a target column and feature columns.

When data is selected, we have established the predictive problem. To start working on the solution, we have to configure model engine (source of model or model architecture), and then configure model hyperparameters. Datagrok automatically builds a UI that updates the model when parameters are changed. So it allows to use Datagrok as interactive playground for modeling.

### Model autoselection and autoconfiguration

Having access to data, Datagrok can suggest models based on the data structure. It suggests a list of models that match the feature columns, but also selects the best matching model from the list.


### Tips as you go

Datagrok analyses data that is used for the prediction and model predictions themselves.
Based on the situation, the platform can suggest data transformations or changes to model configurations or give insights.

Examples of such behavior include:
* If data is imbalanced, Datagrok shows the imbalanced columns
* For data with missing values, Datagrok suggests ignoring rows with missing values or using Missing Values Imputation
* For binary classification, Datagrok analyses false positive and false negative errors.

| Class imbalance | Ignore missing values | Low precision |
|---|---|---|
| ![](./tip-class-imbalance.png) | ![](./tip-impute-missing.png) | ![](./tip-low_accuracy.png) |


### Model comparison

Datagrok allows users to play with models by changing parameters on the go. To avoid losing good parameters found during testing, models can be saved to the model comparison tool.

Model comparison tool saves model parameters, so they can be returned to later. So user can get to the best model configurations and save only it.

## Visualizations

Datagrok supports an extensive list of visualizations, that are shown based on the problem context. For the table, we use $X$ as a notation for the list of columns, and $y$ for the target column. 

| Visualization                | Description                                                                 | Showed for               | Example | Read more                                                                                  |
|------------------------------|-----------------------------------------------------------------------------|--------------------------|---------|-------------------------------------------------------------------------------------------|
| Scatter Plot (regression)    | Visualizes predicted vs actual values.                                      | Regression problems ($y$ is numerical)     | ![](./interactive-scatterplot-numerical.png)  |                                                                                           |
| Residuals plot               | The difference between prediction and the real value.                | Regression problems ($y$ is numerical)     | ![](./interactive-residuals.png)     | [Errors and residuals](https://en.wikipedia.org/wiki/Errors_and_residuals)                |
| Scatter plot (classification)| Visualizes predicted vs actual values using colors for classes.            | Classification problems ($y$ is categorical) | ![](./interactive-scatterplot-categorical.png)      |                                                                                           |
| Roc Curve                    | Trade-offs between true/false positive rates for classifiers.     | Binary classification problems ($y$ is categorical with 2 classes) | ![](./interactive-roc-curve.png)      | [ROC curve](https://en.wikipedia.org/wiki/Receiver_operating_characteristic)              |
| PC Plot                      | Parallel coordinates plot.         | $X$ has more than 2 and less than 11 numerical columns                       | ![](./interactive-pcplot.png)      |                                                                                           |
| Statistics                   | Statistics of the $y$ column                          | Regression problems ($y$ is numerical)          | ![](./interactive-statistics.png)      |                                                                               |
| Distribution      |            Distribution of the predicted values, and distribution of values in $y$.      | Regression problems ($y$ is numerical)                       | ![](./interactive-distributions.png)      |                                                                                  |
| Confusion matrix             | Each row of the matrix represents the instances in an actual class. Each column represents the instances in a predicted class. Also shows relevant metrics.  | Classification problems ($y$ is categorical)  | ![](./interactive-confusion-matrix.png)     | [Confusion matrix](https://en.wikipedia.org/wiki/Confusion_matrix)                        |
| Correlation plot             | Shows feature interdependencies                 | $X$ has more than 10 numerical columns      |    ![](./interactive-corrplot.png) |                                                                                         |
| Performance metrics          | Aggregated performance statistics (e.g., R², accuracy).                    | Always          | ![](./interactive-performance.png)     |                                                                                           |
| Wrong predictions            | Highlights incorrectly classified instances of the data. | Classification problems ($y$ is categorical) | ![](./interactive-mismatches.png)      |                                                                                           |
| Custom visualizations        | Package-defined visualizations for specific models.                    | Specific models          | --      |                                                                                           |

### Custom visualizations

For the custom models, it is possible to define custom viewers using the JS API. It requires a setup by the rules to similar described in [custom models guide](./custom-machine-learning-models.md).

```js

//name: visualize
//meta.mlname: $MODEL_NAME
//meta.mlrole: visualize
//input: dataframe df
//input: column targetColumn
//input: column predictColumn
//input: dynamic model
//output: dynamic widget
export async function visualize(df: DG.DataFrame, targetColumn: DG.Column, predictColumn: DG.Column, model: any): Promise<any> {
  let view : DG.JSViewer = new DG.JSViewer();
  return view.root;
}
```

Such viewers will be automatically added to the interactive training view