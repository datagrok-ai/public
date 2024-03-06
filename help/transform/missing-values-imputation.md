---
title: "Missing values imputation"
---

Use [imputation](https://en.wikipedia.org/wiki/Imputation_\(statistics\)) to replace missing values in a dataframe:

1. On the **Top Menu**, select **ML > Missing Values Imputation...**. A dialog opens.
2. In the dialog, specify the columns with missing values you want to impute (in the `Impute` field) and the columns that should be used for finding neighbors (in the `Using` field). The imputed value is a weighted average of the corresponding values in the specified number of neighbors. You can also select a distance metric, neighbor count, and decide whether to replace missing values or create a new column with imputed results.
3. Click **Run** to execute.

![add-to-workspace](missing-values-imputation.gif)

Datagrok imputes missing values using the k-nearest neighbors method ([k-NN](https://public.datagrok.ai/js/samples/domains/data-science/missing-values-imputation)).

See also:

* [Statistical functions](https://datagrok.ai/help/transform/functions/stats-functions)
* [Recipe Editor](https://datagrok.ai/help/transform/recipe-editor)
