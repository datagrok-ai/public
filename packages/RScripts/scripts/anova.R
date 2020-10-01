#name: ANOVA
#description: One-way ANOVA
#reference: https://en.wikipedia.org/wiki/One-way_analysis_of_variance
#language: r
#tags: demo, analysis, anova
#sample: beer.csv
#input: dataframe table [Input data table]
#input: column categories {type:categorical} [Categories column name]
#input: column variable {type:numerical} [Variable column name]
#output: graphics residualsVsFiltered
#output: graphics normals
#output: graphics scales

result <- aov(table[[variable]] ~ table[[categories]], table)

residualsVsFiltered = plot(result, 1)
normals = plot(result, 2)
scales = plot(result, 3)

print(residualsVsFiltered)
print(normals)
print(scales)
