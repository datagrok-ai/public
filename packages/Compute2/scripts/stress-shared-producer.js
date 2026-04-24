//name: StressSharedProducer
//language: javascript
//tags: stress
//output: dataframe sharedDf {viewer: Grid(allowEdit: false)}

const n = 100;
sharedDf = DG.DataFrame.create(n);
sharedDf.name = 'Shared Data';
sharedDf.columns.addNewInt('id').init((i) => i);
sharedDf.columns.addNewFloat('weight').init(() => Math.random() * 100);
sharedDf.columns.addNewString('label').init((i) => `item_${i}`);
sharedDf.columns.addNewBool('flag').init(() => Math.random() > 0.5);
sharedDf.columns.addNewDateTime('ts').init(() => Date.now() - Math.random() * 1e9);
