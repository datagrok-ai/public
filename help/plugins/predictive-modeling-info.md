<!-- TITLE: Predictive Model Information -->
<!-- SUBTITLE: -->

# Predictive Model Information

Predictive model parameters description.

## Method

Current predictive model method.

### Available classification models

For R Caret via OpenCPU:

| Method    | Model                                                     |
|-----------|-----------------------------------------------------------|
| rf        | Random Forest                                             |
| gbm       | Stochastic Gradient Boosting Machine                      |
| svmLinear | Support Vector Machines with Linear Kernel                |
| svmRadial | Support Vector Machines with Radial Basis Function Kernel |

For H2O:

| Method                       | Model                           |
|------------------------------|---------------------------------|
| Deep Learning                | Deep Learning (Neural Networks) |
| Distributed Random Forest    | DRF                             |
| Generalized Linear Model     | GLM                             |
| Gradient Boosting Machine    | GBM                             |
| Naive Bayes Classifier       |                                 |
| K-Means Clustering           |                                 |
| Principal Component Analysis | PCA                             |


## MSE

The mean squared error (MSE) or mean squared deviation (MSD) of an estimator measures the average of the 
squares of the errors or deviations—that is, the difference between the estimator and what is estimated.

## RMSE

The root-mean-square error (RMSE) is a frequently used measure of the differences between values (sample 
and population values) predicted by a model or an estimator and the values actually observed. The RMSE 
represents the sample standard deviation of the differences between predicted values and observed values.

## NOBS

Number of observations.

## Correlation (Corr or r2)

Correlation is any of a broad class of statistical relationships involving dependence, though in common 
usage it most often refers to the extent to which two variables have a linear relationship with each other. 
The most familiar measure of dependence between two quantities is the Pearson product-moment correlation 
coefficient, or "Pearson's correlation coefficient", commonly called simply "the correlation coefficient". 
It is obtained by dividing the covariance of the two variables by the product of their standard deviations.  

## Log Loss

Logarithmic loss measures the performance of a classification model where the
prediction input is a probability value between 0 and 1. The goal of our machine learning models is to minimize 
this value. A perfect model would have a log loss of 0. Log loss increases as the predicted probability diverges 
from the actual label. So predicting a probability of .012 when the actual observation label is 1 would be
bad and result in a high log loss.

## AUC

The area under the ROC curve (AUC) is equal to the probability that a classifier will rank a randomly chosen
positive instance higher than a randomly chosen negative one (assuming 'positive' ranks higher than 'negative').

## Gini

The Gini coefficient (sometimes expressed as a Gini ratio or a normalized Gini index) measures the inequality 
among values of a frequency distribution. A Gini coefficient of zero expresses perfect equality, 
where all values are the same.

## ROC curve

Receiver operating characteristic curve, i.e. ROC curve, is a graphical plot that illustrates the diagnostic 
ability of a binary classifier system as its discrimination threshold is varied.  
  
The ROC curve is created by plotting the true positive rate (TPR) against the false positive rate (FPR) at 
various threshold settings. The true-positive rate is also known as sensitivity, recall or probability of 
detection in machine learning. The false-positive rate is also known as the fall-out or probability of false 
alarm and can be calculated as (1 - specificity). The ROC curve is thus the sensitivity as a function of fall-out.  
  
ROC analysis provides tools to select possibly optimal models and to discard suboptimal ones independently from 
(and prior to specifying) the cost context or the class distribution. ROC analysis is related in a direct and 
natural way to cost/benefit analysis of diagnostic decision making.  
  
## Confusion matrix

Confusion matrix, also known as an error matrix, is a specific table layout that allows visualization of 
the performance of an algorithm, typically a supervised learning one (in unsupervised learning it is usually 
called a matching matrix). Each column of the matrix represents the instances in a predicted class while each 
row represents the instances in an actual class (or vice versa). The name stems from the fact that it makes it 
easy to see if the system is confusing two classes (i.e. commonly mislabelling one as another).  
  
It is a special kind of contingency table, with two dimensions ("actual" and "predicted"), and identical 
sets of "classes" in both dimensions (each combination of dimension and class is a variable in the contingency table).  
  
See also:

  * [Correlation and dependence](https://en.wikipedia.org/wiki/Correlation_and_dependence)
  * [Root-mean-square error](https://en.wikipedia.org/wiki/Root-mean-square_deviation)
  * [Pearson correlation coefficient](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient)
  * [Receiver operating characteristic (ROC)](https://en.wikipedia.org/wiki/Receiver_operating_characteristic)
  * [Confusion matrix](https://en.wikipedia.org/wiki/Confusion_matrix)
  * [Gini coefficient](https://en.wikipedia.org/wiki/Gini_coefficient)
  * [Log Loss](http://wiki.fast.ai/index.php/Log_Loss)
