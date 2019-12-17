#name: Test
#language: r
#tags: test
#input: int i
#input: double d
#input: bool b
#input: string s
#input: datetime dt
#input: dataframe df
#input: column col
#input: column_list cols
#output: int ri
#output: double rd
#output: bool rb
#output: string rs
#output: datetime rdt
#output: dataframe rdf
#output: dataframe rdfj {action:join(df)}
#test: [12, -30.0, true, "dmc", "1985-10-25 00:00:00", TestData("cars"), "cylinders", ["volume", "price"]]: [6, 30.0, false, "dmcdmc", "1955-11-05 00:00:00", "rdf (6 rows, 3 columns)", "cars (6 rows, 6 columns)"]
ri = i / 2
rd = d + 60
rb = !b
rs = paste(s, s, sep = '')
rdt = as.POSIXlt(dt) - 10947*86400
rdf = df[c(col, cols)]
rdfj = as.data.frame(rowSums(rdf, dims = 1))
colnames(rdfj) = c('metric')
