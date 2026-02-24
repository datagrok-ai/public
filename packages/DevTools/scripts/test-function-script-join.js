//name: TestFunctionScriptJoin
//language: javascript
//meta.vectorFunc: true
//input: dataframe data
//input: column<int> col1
//input: column<string> col2
//input: column<double> col3
//input: list<string> out { optional: true }
//output: dataframe res { action: join(data) }

function createNewCol1() {
  const res1 = DG.Column.int('joinedCol1', col1.length);
  res1.init((i) => col1.getNumber(i) + 1);
  return res1;
}
function createNewCol2() {
  const res2 = DG.Column.string('joinedCol2', col2.length);
  res2.init((i) => col2.get(i) + ' joined');
  return res2;
}
function createNewCol3() {
  const res3 = DG.Column.float('joinedCol3', col3.length);
  res3.init((i) => col3.getNumber(i) + 10.5);
  return res3;
}
const colCreationFuncs = {
  'joinedCol1': createNewCol1,
  'joinedCol2': createNewCol2,
  'joinedCol3': createNewCol3,
};

let colList = [];
if (out == undefined || out.length === 0)
  colList = [createNewCol1(), createNewCol2(), createNewCol3()];
else {
  for (const colName of out) {
    if (colCreationFuncs[colName] != undefined)
      colList.push(colCreationFuncs[colName]());
  }
}
res = colList.length > 0 ? DG.DataFrame.fromColumns(colList) : DG.DataFrame.create(col1.length);
