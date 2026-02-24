//name: TestFunctionScriptJoinColumnList
//language: javascript
//meta.vectorFunc: true
//input: dataframe data
//input: column_list columns
//input: list<string> out { optional: true }
//output: dataframe res { action: join(data) }

let cols = columns.toList();
let colCreationFuncs = {};
for (let i = 0; i < cols.length; i++) {
  const col = cols[i];
  const idx = i;
  colCreationFuncs[`joinedCol${idx + 1}`] = () => {
    if (col.type === DG.COLUMN_TYPE.INT || col.type === DG.COLUMN_TYPE.FLOAT) {
      const r = DG.Column.float(`joinedCol${idx + 1}`, col.length);
      r.init((j) => col.getNumber(j) + idx + 1);
      return r;
    } else {
      const r = DG.Column.string(`joinedCol${idx + 1}`, col.length);
      r.init((j) => col.get(j) + ` joined${idx + 1}`);
      return r;
    }
  };
}

let colList = [];
if (out == undefined || out.length === 0) {
  for (const name of Object.keys(colCreationFuncs))
    colList.push(colCreationFuncs[name]());
} else {
  for (const colName of out) {
    if (colCreationFuncs[colName] != undefined)
      colList.push(colCreationFuncs[colName]());
  }
}
res = colList.length > 0 ? DG.DataFrame.fromColumns(colList) : DG.DataFrame.create(data.rowCount);
