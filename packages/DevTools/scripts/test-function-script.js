//name: TestFunctionScript
//language: javascript
//meta.vectorFunc: true
//input: column<int> col1
//input: column<string> col2
//input: column<double> col3
//input: list<string> out { optional: true }
//output: dataframe res

function createNewCol1() {
  const res1 = DG.Column.int('newCol1', col1.length);
  res1.init((i) => col1.getNumber(i) + 1);
  return res1;
}
function createNewCol2() {
  const res2 = DG.Column.string('newCol2', col2.length);
  res2.init((i) => col2.get(i) + ' and 123');
  return res2;
}
function createNewCol3() {
  const res3 = DG.Column.float('newCol3', col3.length);
  res3.init((i) => col3.getNumber(i) + 5.5);
  return res3;
}
const colCreationFuncs = {
  'newCol1': createNewCol1,
  'newCol2': createNewCol2,
  'newCol3': createNewCol3,
};

console.log(out);
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
