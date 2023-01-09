#name: Scatter Plot Julia
#language: julia
#tags: demo, viewers, hide-suggestions
#input: dataframe t
#input: column xColumnName
#input: column yColumnName
#input: column colorColumnName
#output: graphics

using Plots;
colors = [RGB(x, x, x) for (x) in t[colorColumnName]]
scat = scatter(t[xColumnName], t[yColumnName], markercolor=colors, xlabel=xColumnName, ylabel=yColumnName, legend=false)
display(scat)
