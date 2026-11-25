//name: ViewersHookModel
//language: javascript
//description: RichFunctionView viewers hook test
//input: int n = 3
//output: dataframe df {viewer: Grid(allowEdit: false)}
//editor: Compute2:RichFunctionViewEditor
//meta.viewersHook: Compute2:SlowViewersHookMaker
//tags: test

df = DG.DataFrame.create(n);
