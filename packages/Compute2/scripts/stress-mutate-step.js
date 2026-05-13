//name: StressMutateStep
//language: javascript
//tags: stress
//input: dataframe df {viewer: Grid(allowEdit: false)}
//input: dataframe sharedDf {viewer: Grid(allowEdit: false)}
//output: dataframe dfOut {viewer: Grid(allowEdit: false)}

const n = df ? df.rowCount : 100;
const cols = [];

cols.push(DG.Column.fromType('int', 'id', n).init((i) => df ? df.get('id', i) : i));
cols.push(DG.Column.fromType('double', 'val1', n).init((i) => df ? df.get('val1', i) : Math.random()));
cols.push(DG.Column.fromType('double', 'val2', n).init((i) => df ? df.get('val2', i) : Math.random()));
cols.push(DG.Column.fromType('string', 'cat1', n).init((i) => df ? df.get('cat1', i) : `cat_${i % 10}`));
cols.push(DG.Column.fromType('string', 'cat2', n).init((i) => df ? df.get('cat2', i) : `grp_${i % 5}`));
cols.push(DG.Column.fromType('bool', 'active', n).init((i) => df ? df.get('active', i) : i % 2 === 0));
cols.push(DG.Column.fromType('int', 'count', n).init((i) => df ? df.get('count', i) : Math.floor(Math.random() * 1000)));
cols.push(DG.Column.fromType('double', 'score', n).init((i) => df ? df.get('score', i) : Math.random() * 100));
cols.push(DG.Column.fromType('datetime', 'timestamp', n).init(() => Date.now() - Math.random() * 1e9));
cols.push(DG.Column.fromType('double', 'noise', n).init(() => Math.random()));

dfOut = DG.DataFrame.fromColumns(cols);
dfOut.name = 'Mutated';

// Randomly mutate 5 values
for (let m = 0; m < 5; m++) {
  const row = Math.floor(Math.random() * n);
  const colIdx = Math.floor(Math.random() * 3) + 1; // mutate val1, val2, or score
  const col = dfOut.columns.byIndex(colIdx);
  if (col.type === 'double')
    col.set(row, Math.random() * 100);
  else if (col.type === 'int')
    col.set(row, Math.floor(Math.random() * 1000));
}
