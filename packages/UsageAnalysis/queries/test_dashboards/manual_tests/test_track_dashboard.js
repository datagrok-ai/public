//language: javascript
//input: dataframe result
//output: dataframe out

async function postprocess() {
  let nameCategories = result.col('batch_name').categories;
  let batchNames = [];
  
  for (var i = 0; i < result.col('batch_name').length && nameCategories.length > 0; i++) {
    let nm =result.col('batch_name').get(i);
    if (!batchNames.includes(nm))
    {
      nameCategories = nameCategories.filter(e => e !== nameCategories.includes(nm));
      batchNames.push(nm)
    }
  }
  let order = ['test', ...(batchNames.map(e=> `${e} concat unique(status)`))];
  let pivot = result
    .groupBy(['test'])
    .pivot('batch_name')
    //.avg('date_time')
    .add(DG.STR_AGG.CONCAT_UNIQUE, 'status')
    .add(DG.STR_AGG.CONCAT_UNIQUE, 'result')
    .aggregate();

  function replaceColumn(prefix, type, newType) {
    let col = pivot.columns.byName(prefix + type);
    col?.setTag('friendlyName', prefix + newType);
  }

  pivot.columns.byName('test').semType = 'autotest'; 
  for (let i = 0; i < batchNames.length; i++) {
    let batchName = batchNames[i]; 
    replaceColumn(batchNames[i], ' concat unique(status)', '');
    replaceColumn(batchNames[i], ' concat unique(result)', ' verdict');
    let statusCol = pivot.col(`${batchNames[i]} concat unique(status)`);
    if (statusCol)
      statusCol.meta.colors.setCategorical({
        'passed': '#2ca02c',
        'skipped': '#ffa500',
        'failed': '#8e342a',
      })
    let col = pivot.col(`${batchName} concat unique(result)`);
    for (var j = 0; j < col.length; j++) {
      let tickets = col.getString(j);
      col.set(j, tickets.replaceAll('\n', ','));
    }
    order.push(`${batchNames[i]} concat unique(result)`);
  }
 pivot.columns.setOrder(order)

  pivot.name = 'test track dashboard';
  result.rows.filter(x => x.batch_name === result.rows.get(0).batch_name)
  out = pivot;
}

if (result.rowCount == 0) {
  out = DG.DataFrame.fromColumns([
    DG.Column.fromType('string', 'test', 1),
    DG.Column.fromType('string', 'owner', 1)
  ]);
}
else
  await postprocess();