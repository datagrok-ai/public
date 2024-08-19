//name: JavaScript Graphics
//description: graphics output column input
//language: javascript
//input: dataframe df
//input: column xName {type:numerical; allowNulls:false} [Column x]
//input: column yName {type:numerical; allowNulls:false} [Column y]
//output: graphics result

view = grok.shell.addTableView(df);
result = df.plot.grid({x: xName, y: yName});
