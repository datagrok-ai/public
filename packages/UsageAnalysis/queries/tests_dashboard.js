//language: javascript
//input: dataframe result
//output: dataframe out  

let builds = result.col('build_index').categories;
let pivot = result
  .groupBy(['test', 'owner'])
  .pivot('build_index')
  .avg('duration')
  .avg('build_date')
  .add(DG.STR_AGG.CONCAT_UNIQUE, 'status')
  .add(DG.STR_AGG.CONCAT_UNIQUE, 'result')
  .aggregate(); 

function generateCommonValueFormula(operation, value) {
  var result = '';
  if (builds.length == 1)
    return '${' + builds[0] + ' concat unique(status)}' + value;
  for (var i = 0; i + 2 < builds.length; i++)
    result += operation + '(${' + builds[i] + ' concat unique(status)}' + value + ',';
  result += operation + '(${' + builds[builds.length - 2] + ' concat unique(status)}' + value + ',${' +  builds[builds.length - 2] + ' concat unique(status)}' + value + ')';
  for (var i = 0; i + 2 < builds.length; i++)
    result += ')';
  return result;
}

pivot.columns.addNewCalculated('stable', generateCommonValueFormula("And", '== "passed"'));
pivot.columns.addNewCalculated('failing', generateCommonValueFormula("Or", '== "failed"'));
pivot.columns.addNewCalculated('flapping', "And(${failing}, " + generateCommonValueFormula("Or", '== "passed"') + ")");


out = pivot;
