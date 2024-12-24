// Create dataframe by adding and initializing columns
//tags: DataFrame, Column, construction

let t = DG.DataFrame.create(10)
t.columns.addNewInt('int').init((i) => i);
t.columns.addNewFloat('float').init((i) => i / 10);
t.columns.addNewQnum('qnum').init((i) => DG.Qnum.create(i / 10, DG.QNUM_LESS));
t.columns.addNewString('string').init((i) => 'str ' + i);
t.columns.addNewBool('bool').init((i) => i % 2 === 0);
t.columns.addNewDateTime('datetime').init((i) => dayjs());
t.columns.addNewBytes('bytes').init((i) => new Uint8Array(42));

grok.shell.addTableView(t);
