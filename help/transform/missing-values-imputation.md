<!-- TITLE: Missing Values Imputation -->
<!-- SUBTITLE: -->

# Missing Values Imputation

Imputation is the process of replacing missing data with substituted values.  
  
## Algorithm

Missing values imputation algorithm is based on **k-nearest neighbors algorithm (k-NN)**. This
algorithm is a non-parametric method used for classification and regression. Both for classification
and regression, it can be useful to assign weight to the contributions of the neighbors, so that the
nearer neighbors contribute more to the average than the more distant ones. For example, a common
weighting scheme consists in giving each neighbor a weight of **1/d**, where **d** is the distance
to the neighbor.

## Run

  * Open table
  * Run from menu: **Tools** | **Data Science** | **Missing Values Imputation...**
  * Select source table
  * Select columns that contains missing values **"impute"**
  * Select data rows **"data"**
  * Set **"Number of nearest neighbours"**
  * Run missing values imputation. Result will replace all missing values in all columns and rows.  
  
See also:
  * [Imputation](https://en.wikipedia.org/wiki/Imputation_\(statistics\))
  * [k-nearest neighbors algorithm](https://en.wikipedia.org/wiki/K-nearest_neighbors_algorithm)
  
Sample: 
  * [Missing Values Imputation](https://public.datagrok.ai/js/samples/domains/data-science/missing-values-imputation)
