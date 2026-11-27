# Machine learning and statistics

Playlist: **Machine learning and statistics**. Feature name of every clip: `How-to · ML`. Features in
`public/packages/Tutorials/bdd/features/howto/10-ml-and-statistics/`. Videos in
`~/Desktop/grok/Core/howto-videos/10-ml-and-statistics/`.

## Plan

One idea per clip, for a scientist who knows what PCA or a regression is but not where Datagrok keeps it. Only
tools that run in the browser (the EDA package's TypeScript and Crux/WASM, Forge's in-browser methods, Dendrogram's
workers); Python scripts, notebooks and containers stay out. Cluster Data (k-means) and Impute Missing Values are in
the "Work with data" playlist; chemical space and activity cliffs in the molecules one; the scatter plot's
regression line and the correlation plot in the viewer playlists.

1. What the ML menu does: a walk down its groups (demog).
2. PCA: many numeric columns into two components, plotted and colored by a category (winequality).
3. UMAP / t-SNE: a 2D map of the rows from Reduce Dimensionality, colored by a category (winequality).
4. Train a model and apply it: ML > Models > Train Model (Forge), target, features, method, the scores on held-out
   rows, another method, Save, where saved models are; then ML > Models > Apply Model on a table without the target
   column (iris; nine iris rows). Merged on 2026-11-22 from "Train a model" and "Apply a model to new data" (0:33).
6. Multivariate analysis (PLS): what drives price, the charts it docks (cars).
7. ANOVA: does a number differ between groups, and how to read the answer (demog-1000).
8. Pareto front: the best trade-offs between price and horsepower, and what green and grey mean on three cars
   (cars). Extended on 2026-11-22 from 0:32.
9. Hierarchical clustering: a tree beside the grid, cut into clusters (iris).

Left out: the t-test and Control Comparisons (the same dialog and result as ANOVA, for two groups or against a
control), DBSCAN and MCL (other clustering methods; k-means is in "Work with data"), Random Data, Notebooks
(Jupyter).

## 01 - What the ML menu does (0:42)

Title: The ML menu in Datagrok: a quick tour

Description: A walk down the ML menu of a table view, one group at a time: Analyze (PCA, PLS, multivariate
analysis and group comparisons such as the t-test and ANOVA), Reduce Dimensionality (UMAP and t-SNE), Impute Missing
Values, Pareto Front, Cluster (k-means, DBSCAN, MCL, hierarchical clustering) and Models (train a predictive model,
apply a saved one). Nothing is run here; most of these tools have a video of their own in this playlist. Uses the
demog demo dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/10-ml-and-statistics/01-the-ml-menu.feature`

## 02 - Reduce many columns with PCA (1:08)

Title: PCA in Datagrok: many columns into two

Description: Principal component analysis turns many numeric columns into a few new ones that keep most of the
differences between rows. ML > Analyze > PCA opens a dialog where you pick the feature columns and the number of
components; OK adds PC1 and PC2 to the table. A scatter plot of PC1 against PC2, colored by wine type, shows red
and white wines on opposite sides, though PCA never saw the type. Uses the winequality demo dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/10-ml-and-statistics/02-pca.feature`

## 03 - Map similar rows with UMAP (0:59)

Title: UMAP and t-SNE in Datagrok: a map of similar rows

Description: ML > Reduce Dimensionality places every row on a 2D map so that rows with similar values end up close
together. You pick the columns to compare and the method, UMAP or t-SNE; OK computes the map in the browser, adds
the coordinates to the table and plots them, colored by the clusters it finds. Colored by wine type instead, the red
and white wines land in different parts of the map. The map's shape changes from run to run; which rows end up
near each other does not. Uses the winequality demo dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/10-ml-and-statistics/03-map-rows-with-umap.feature`

## 04 - Train a model and apply it (1:34)

Title: Train a predictive model and apply it to new data in Datagrok

Description: ML > Models > Train Model opens the training view: the column to predict, the columns it learns
from, and the method. It trains right away in the browser and shows its scores on rows held out by 5-fold
cross-validation; picking another method retrains it, and Runs keeps every training for comparison. Save stores
the model under a name, and Browse > Platform > Machine Learning lists your saved models. Then, in a table of new
measurements, ML > Models > Apply Model picks the saved model, matches its inputs to the columns of the same
name, and adds the prediction as a new column. Uses the iris demo dataset and nine new flowers without a species.

Feature: `public/packages/Tutorials/bdd/features/howto/10-ml-and-statistics/04-train-a-model.feature`

## 05 - Find what drives a value with multivariate analysis (1:12)

Title: Multivariate analysis (PLS) in Datagrok: what drives a value

Description: ML > Analyze > Multivariate Analysis fits a partial least squares model of one column from many others.
You choose the column to predict, the columns that might explain it, the number of components and a column to label
the points; RUN adds the predicted values and docks the charts that explain the model: observed against predicted,
variable importance, loadings, scores and explained variance. Uses the cars demo dataset, where engine size and
horsepower matter most for price.

Feature: `public/packages/Tutorials/bdd/features/howto/10-ml-and-statistics/05-multivariate-analysis.feature`

## 06 - Compare groups with ANOVA (0:42)

Title: One-way ANOVA in Datagrok: compare a value across groups

Description: ML > Analyze > Group Comparison > ANOVA tests whether the average of a numeric column differs between
the groups of a category. You pick the category, the number to compare and the significance level; Run docks a box plot with the
conclusion and its p-value on top, and a table with the F statistic, its critical value and the p-value. Uses the
demog-1000 demo dataset: age differs between the six disease groups.

Feature: `public/packages/Tutorials/bdd/features/howto/10-ml-and-statistics/06-anova.feature`

## 07 - Find the best trade-offs with a Pareto front (0:48)

Title: Pareto front in Datagrok: the best trade-offs

Description: When goals pull against each other, the Pareto front marks the rows that no other row beats on all of
them at once. ML > Pareto Front opens an Optimize panel where you set each column to minimize or maximize, a plot of
the front and a parallel coordinates plot. With price minimized and horsepower maximized on the cars demo dataset,
the green cars are the ones no other car beats on both: the Mercedes is grey because the Alfa Romeo has more power
for about half the price, while the Porsche, the dearest car, is green because none has more power.

Feature: `public/packages/Tutorials/bdd/features/howto/10-ml-and-statistics/07-pareto-front.feature`

## 08 - Cluster rows into a tree (1:05)

Title: Hierarchical clustering in Datagrok: a tree beside the grid

Description: ML > Cluster > Hierarchical Clustering joins similar rows step by step and draws the result as a tree
next to the grid, one leaf per row. You pick the feature columns, the distance and the linkage; right-click the tree
and Assign Clusters cuts it into the number of groups you choose and adds a Cluster column, and Color Coding >
Categorical in the column header's menu colors the groups. Uses the iris demo dataset, with the rows made short so
all 150 flowers fit: setosa gets a cluster of its own, while versicolor and virginica overlap.

Feature: `public/packages/Tutorials/bdd/features/howto/10-ml-and-statistics/08-hierarchical-clustering.feature`
