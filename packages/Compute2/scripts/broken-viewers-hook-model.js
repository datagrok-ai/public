//name: BrokenViewersHookModel
//language: javascript
//description: RichFunctionView whose viewers hook can't be made
//input: int n = 2
//output: dataframe df {viewer: Grid()}
//editor: Compute2:RichFunctionViewEditor
//meta.viewersHook: Compute2:FailingViewersHookMaker
//tags: test

df = DG.DataFrame.create(n);
