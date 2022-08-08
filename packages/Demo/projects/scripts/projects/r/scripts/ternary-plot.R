#name: Ternary Plot
#language: r
#tags: demo, viewers
#sample: feldspar.csv
#input: dataframe t
#input: column topColumnName {type:numerical;auto}
#input: column leftColumnName {type:numerical;auto}
#input: column rightColumnName {type:numerical;auto}
#input: int pointSize = 3
#output: graphics

require(ggtern)

data <- data.frame(top=t[[topColumnName]], left=t[[leftColumnName]], right=t[[rightColumnName]])
plotTernary = ggtern(data=data,aes(top,left,right)) +
  Tlab(topColumnName) + Llab(leftColumnName) + Rlab(rightColumnName) +
  stat_density_tern(
    geom='polygon',
    aes(fill=..level..),
    bins=5,
    color='grey') +
  geom_point(size=pointSize)

print(plotTernary)
