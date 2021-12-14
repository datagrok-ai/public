#name: t-test Julia
#description: Welch's t-test
#reference: https://en.wikipedia.org/wiki/Welch%27s_t-test
#language: julia
#tags: demo, hide-suggestions
#sample: TSLA.csv
#input: dataframe data [Input data table]
#input: column x {type:numerical, allowNulls:false} [X axis column name]
#input: column y {type:numerical, allowNulls:false} [Y axis column name]
#output: double pValue [P-value of t-statistics]

using HypothesisTests

pValue = pvalue(UnequalVarianceTTest(disallowmissing(data[x]), disallowmissing(data[y])));
