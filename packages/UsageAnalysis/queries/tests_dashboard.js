//language: javascript
//input: dataframe result
//output: dataframe out  

let builds = result.col('build').categories;
let pivot = result
  .groupBy(['test'])
  .pivot('build')
  .avg('duration')
  .add(DG.STR_AGG.CONCAT_UNIQUE, 'status')
  .add(DG.STR_AGG.CONCAT_UNIQUE, 'result')
  .aggregate(); 

function generateCommonValueFormula(operation, value) {
  var result = '';
  if (builds.length == 1)
    return builds[0] + value;
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
