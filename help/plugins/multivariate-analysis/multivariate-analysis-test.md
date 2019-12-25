<!-- TITLE: Tests: Multivariate Analysis -->
<!-- SUBTITLE: -->

# Tests: Multivariate Analysis

[Multivariate analysis (MVA)](../plugins/multivariate-analysis/pls.md) is based on the statistical principle of multivariate statistics.


## Testing scenarios

1. Open table "demog"

1. Open "Partial Least Squares (PLS)" from **Tools | Data Science | Multivariate Analysis (PLS)...**

1. Set value of the field "Outcome" as a column which have empty rows (HEIGHT)
   * After selecting a column with nulls for "Outcome" field, a warning is displayed 
   * "OK" does not available for clicking

1. Select columns with nulls for "Features". (WEIGHT)
   * After selecting a column with nulls for, a warning is displayed 
   * "OK" does not avaible for clicking
   * In warning there is an opportunity to use ["MissingValuesImputation"](../dialogs/missing-values-imputation.md)

1. Click on ["MissingValuesImputation"](../dialogs/missing-values-imputation.md) on the warning about nulls values.
Use ["MissingValuesImputation"](../dialogs/missing-values-imputation.md) instrument for input missing values to WEIGNT column.
   * After clicking on the ["MissingValuesImputation"](../dialogs/missing-values-imputation.md), an appropriate dialog was opened
   * After imputation missing values, WEIGHT column became available for selection and does not cause notifications

1. Select another columns with nulls for "Features". (HEIGHT)
   * After selecting a column with nulls for, a warning is displayed. 
   * "OK" does not avaible for clicking

1. Without closing "PLS" dialog, open **Select | "Missing Values"**, then select the "HEIGHT" column and delete the selected empty rows.
   * After deleting rows with missing values, HEIGHT column became available for selection and does not cause notifications

1. Check "Features" and "Outcome" fields for using only numeric columns.
   * Only numeric columns available for select for "Features" and "Outcome" fields.

1. Select valid columns for "Features" field (HEIGHT, WEIGHT after deleting nulls) and for "Outcome" field (AGE column without nulls). 
   Run the "[PLS](../plugins/multivariate-analysis/pls.md)".

See also: 
  * [PLS](../plugins/multivariate-analysis/pls.md)
  * [Correlation Loadings plot](../plugins/multivariate-analysis/plots/correlation-loadings.md)
  * [Explained Variance plot](../plugins/multivariate-analysis/plots/explained-variance.md)
  * [Regression Coefficients plot](../plugins/multivariate-analysis/plots/regression-coefficients.md)
  * [Predicted vs. Reference plot](../plugins/multivariate-analysis/plots/predicted-vs-reference.md)
  * [Multivariate analysis Auto Test](multivariate-analysis-test.side)
