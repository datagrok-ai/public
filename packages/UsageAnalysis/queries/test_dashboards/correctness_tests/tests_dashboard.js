//language: javascript
//input: dataframe result
//output: dataframe out

function getFirstOfCol(col) {
  if (col)
    return col.categories.filter(x => x != '')[0];
  return '';
}


async function postprocess() {

  let builds = result.col('build_index').categories.map((x) => parseInt(x));
  let buildNames = result.col('build').categories.sort();
  let pivot = await result
    .groupBy(['test', 'owner'])
    .pivot('build_index')
    .avg('duration')
    .avg('build_date')
    .first('flaking')
    .add(DG.STR_AGG.CONCAT_UNIQUE, 'status')
    .add(DG.STR_AGG.CONCAT_UNIQUE, 'instance')
    .add(DG.STR_AGG.CONCAT_UNIQUE, 'build_commit')
    .add(DG.STR_AGG.CONCAT_UNIQUE, 'build')
    .add(DG.STR_AGG.CONCAT_UNIQUE, 'result')
    .aggregate();

  function generateCommonValueFormula(operation, suffix, value) {
    let res = '';
    if (builds.length == 1)
      return '${' + builds[0] + suffix + value;
    for (let i = 0; i + 2 < builds.length; i++)
      res += operation + '(${' + builds[i] + suffix + value + ',';
    res += operation + '(${' + builds[builds.length - 2] + suffix + value + ',${' +  builds[builds.length - 1] + suffix + value + ')';
    for (let i = 0; i + 2 < builds.length; i++)
      res += ')';
    return res;
  }
  await pivot.columns.addNewCalculated('stable', generateCommonValueFormula("And", ' concat unique(status)}', '== "passed"'), 'bool');
  await pivot.columns.addNewCalculated('failing', generateCommonValueFormula("Or", ' concat unique(status)}', '== "failed"'), 'bool');
  await pivot.columns.addNewCalculated('flaking', generateCommonValueFormula("Or", ' first(flaking)}', '== true'), 'bool');
  let schemas = await grok.dapi.stickyMeta.getSchemas();
  let schema = schemas.filter((schema) => schema.name == 'Autotests').at(0);
  let meta = await grok.dapi.stickyMeta.getAllValues(schema, pivot.columns.byName('test'));
  pivot.columns.add(meta.col('ignore?'));
  pivot.columns.add(meta.col('ignoreReason'));
  pivot.columns.add(meta.col('lastResolved'));  
  let needs_attention = await pivot.columns.addNewCalculated('needs_attention', "And(${failing}, Not(${ignore?}))", 'bool');

  function replaceColumn(prefix, type, buildName, newType) {
    let col = pivot.columns.byName(prefix + type);
    if (col != null)
      col.name = buildName + newType;
    col?.setTag('friendlyName', buildName + newType);
  }
    
    
  function addTagsToColumn(col, tags) {
    if (col)
      for (let i of Object.keys(tags)) {
          if (tags[i])
            col?.setTag(i, tags[i]);
    
      }
  }

  pivot.columns.byName('owner').semType = 'User';
  pivot.columns.byName('owner').setTag('cell.renderer', 'User');
  pivot.columns.byName('test').semType = 'autotest';
  
  
  let ticketsStatusCol = await pivot.columns.addNewCalculated(`${builds[0]} tickets status`, `UsageAnalysis:getTicketsVerdict(\${${builds[0]} concat unique(result)})`, DG.TYPE.STRING);
  if (ticketsStatusCol)
    ticketsStatusCol.colors.setCategorical({
      'Fixed': '#2ca02c',
      'Partially Fixed (Lowest)': '#ffe51c',
      'Partially Fixed (Low)': '#ffa500',
      'Partially Fixed (Medium)': '#ff5a00',
      'Partially Fixed (Blocker)': '#7b2d24',
      'Wasn\'t Fixed (Lowest)': '#ffe51c',
      'Wasn\'t Fixed (Low)': '#ffa500',
      'Wasn\'t Fixed (Medium)': '#ff5a00',
      'Wasn\'t Fixed (Blocker)': '#7b2d24',
    })
    

  for (let i = 1; i <= builds.length; i++) {
    let buildName = builds[i - 1];
    let data = {build : (getFirstOfCol(pivot.columns.byName(buildName + ' concat unique(build)'))),
     build_commit : (getFirstOfCol(pivot.columns.byName(buildName + ' concat unique(build_commit)'))),
     build_date : (getFirstOfCol(pivot.columns.byName(buildName + ' avg(build_date)'))),
     instance : (getFirstOfCol(pivot.columns.byName(buildName + ' concat unique(instance)')))}
    let colResult = pivot.columns.byName(i + ' concat unique(result)');
    if (colResult !== null)
      colResult.semType = 'stackTrace';
    addTagsToColumn(pivot.getCol(buildName + ' avg(duration)'), data);
    addTagsToColumn(pivot.getCol(buildName + ' concat unique(status)'), data);
    replaceColumn(buildName, ' avg(duration)', i, ' duration');
    replaceColumn(buildName, ' first(flaking)', i, ' flaking');
    replaceColumn(buildName, ' concat unique(result)', i, ' result');
    replaceColumn(buildName, ' concat unique(status)', i, ''); 
  }

  for (let i = 0; i < pivot.rowCount; i++) {
  }

  for (let i = 1; i <= builds.length; i++) {
    let buildName = builds[i - 1];
    pivot.columns.remove(buildName + ' concat unique(build)');
    pivot.columns.remove(buildName + ' concat unique(build_commit)');
    pivot.columns.remove(buildName + ' avg(build_date)');
    pivot.columns.remove(buildName + ' concat unique(instance)');
  }

  pivot.name = '0. Tests Dashboard'; 

 return pivot;
}
let out = null;
if (result.rowCount == 0)  {
  out = DG.DataFrame.fromColumns([
    DG.Column.fromType('string', 'test', 1),
    DG.Column.fromType('string', 'owner', 1)
  ]);
}
else
    out = await postprocess();