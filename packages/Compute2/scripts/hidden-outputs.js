//name: HiddenOutputsTest
//language: javascript
//description: RichFunctionView hidden dataframe output test
//input: double a
//output: dataframe shown
//output: dataframe secret {hidden: true}
//editor: Compute2:RichFunctionViewEditor
//tags: test

shown = DG.DataFrame.fromCsv(`x\n${a}`);
secret = DG.DataFrame.fromCsv(`x\n${a}`);
