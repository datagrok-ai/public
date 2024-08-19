#name: R Column List
#description: column list input
#language: r
#input: dataframe df
#input: column_list cols
#output: dataframe result

result <- df[,names(df) %in% cols]
