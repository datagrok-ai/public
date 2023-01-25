#name: t-test Python
#description: Welch's t-test
#reference: https://en.wikipedia.org/wiki/Welch%27s_t-test
#language: python
#tags: demo, hide-suggestions
#sample: TSLA.csv
#input: dataframe data [Input data table]
#input: column x {type:numerical} [X axis column name]
#input: column y {type:numerical} [Y axis column name]
#output: double pValue [P-value of t-statistics]

from scipy import stats

_, p = stats.ttest_ind(data[[x]], data[[y]], equal_var=False)
pValue = p[0]
