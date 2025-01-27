//language: javascript
//input: dataframe result
//output: dataframe out

async function postprocess() {

  let builds = result.col('build_index').categories;
  let buildNames = result.col('build').categories.sort().reverse();

  // console.log(result.name);
  // console.log(buildNames[0]);

  let pivot = result
    .groupBy(['test','owner'])
    .pivot('build_index')
    .pivot('worker')
    .pivot('browser')
    .avg('build_date')
    .add(DG.STR_AGG.CONCAT_UNIQUE, 'status')
    .aggregate();
  // console.log(pivot.name);


  function generateCommonValueFormula(operation, suffix, value) {
    var result = '';
    if (builds.length == 1)
      return '${' + builds[0] + suffix + value;
    for (var i = 0; i + 2 < builds.length; i++)
      result += operation + '(${' + builds[i] + suffix + value + ',';
    result += operation + '(${' + builds[builds.length - 2] + suffix + value + ',${' +  builds[builds.length - 1] + suffix + value + ')';
    for (var i = 0; i + 2 < builds.length; i++)
      result += ')';
    return result;
  }

  pivot.columns.byName('owner').semType = 'User';
  pivot.columns.byName('owner').setTag('cell.renderer', 'User');
  pivot.columns.byName('test').semType = 'test';

  // await pivot.columns.addNewCalculated('failing', generateCommonValueFormula("Or", ' concat unique(status)}', '== "failed"'), 'bool');
  // await pivot.columns.addNewCalculated('flaking', generateCommonValueFormula("Or", ' first(flaking)}', '== true'), 'bool');
  // await pivot.columns.addNewCalculated('needs_attention', "Or(${flaking},${failing})", 'bool');

  function replaceColumn(prefix, type, buildName, newType) {
    var col = pivot.columns.byName(prefix + type);
    if (col != null)
      col.name = '{' + col?.name;
    col?.setTag('friendlyName', buildName + newType);
  }

  var groups = {};
  for (var i = 1; i <= builds.length; i++) {
    var columns = [];
    var dateColumn = null;
    for (var j = 0; j < 5; j++) {
      for (var l = 0; l < 5; l++) {
        var buildName = buildNames[i - 1];
        replaceColumn(i, ' ' + j + ' ' + l + ' concat unique(status)', buildName + ':' + j, '');
        if (dateColumn == null) {
          dateColumn = pivot.columns.byName(i +  ' ' + j + ' ' + l + ' avg(build_date)');
        } else {
          var col = pivot.columns.byName(i +  ' ' + j + ' ' + l + ' avg(build_date)');
          if (col != null)
            pivot.columns.remove(i +  ' ' + j + ' ' + l + ' avg(build_date)');
        }
        columns.push('{' + i + ' ' + j + ' ' + l + ' concat unique(status)');
        var col = pivot.columns.byName('{' + i + ' ' + j + ' ' + l + ' concat unique(status)');
        if (col !== null) {
          for (var k = 0; k < col.length; k++) {
            if (col.get(k) == 'failed') {
              var failedCol = pivot.columns.getOrCreate(i + ': failed', DG.TYPE.INT);
              var old = failedCol.isNone(k) ? 0 : failedCol.get(k);
              failedCol.set(k, old + 1);
            }
            if (col.get(k) == 'passed') {
              var successCol = pivot.columns.getOrCreate(i + ': success', DG.TYPE.INT);
              var old = successCol.isNone(k) ? 0 : successCol.get(k);
              successCol.set(k, old + 1);
            }
          }
        }
      }
      if (dateColumn != null) {
        dateColumn.name = i + ': date'; 
        dateColumn.setTag('friendlyName', buildName + ': date');
      }
    }
    groups[buildNames[i - 1]] = {
      columns: columns
    };
  }

  pivot.meta.setGroups(groups);
  var sortedNames = pivot.columns.names().sort();
  pivot.columns.setOrder(sortedNames);
    
  // const jsonColumn = pivot.columns.addNewString('duration');
  // jsonColumn.semType = FIT_SEM_TYPE;
  // for (let i = 0; i < pivot.rowCount; i++) {
  //   var chartData = {
  //     series: [{
  //       fitFunction: "linear",
  //       points: [
  //         {x: 1, y: pivot.columns.byName('1 avg(duration)').get(i)},
  //         {x: 2, y: pivot.columns.byName('2 avg(duration)').get(i)},
  //         {x: 3, y: pivot.columns.byName('3 avg(duration)').get(i)},
  //         {x: 4, y: pivot.columns.byName('4 avg(duration)').get(i)},
  //         {x: 5, y: pivot.columns.byName('5 avg(duration)').get(i)},
  //       ]
  //     }]
  //   };
  //   jsonColumn.set(i, JSON.stringify(chartData));
  // }


  pivot.name = 'stress test dashboard';

  out = pivot;

}

if (result.rowCount == 0) 
  out = DG.DataFrame.fromColumns([
    DG.Column.fromType('string', 'test'),
    DG.Column.fromType('string', 'test')
  ]);
else
  await postprocess();