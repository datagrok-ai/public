#name: ListPackages
#description: Body Mass Index
#language: r
#output: dataframe x

x = as.data.frame(installed.packages()[,c(1,3:4)])