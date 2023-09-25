//name: Data stats
//description: Uses Datagrok's built-in viewer to show dataframe stats
//language: javascript
//tags: demo
//input: dataframe inputDf {caption: Input dataframe; viewer: Grid(); category: Input data}
//output: dataframe outputDf {caption: Output dataframe; viewer: Line chart(block: 50) | Scatter plot(block: 50) | Statistics(block: 100); category: Stats}
//editor: Compute:RichFunctionViewEditor
//meta.runOnOpen: true
//meta.runOnInput: true


outputDf = inputDf.clone();
