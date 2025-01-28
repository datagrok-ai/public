//language: javascript
//input: dataframe result
//output: dataframe out


async function postprocess() {

  let builds = result.col('build_index').categories;
  let buildNames = result.col('build').categories.sort().reverse();

  let pivot = result
    .groupBy(['test', 'owner'])
    .pivot('build_index')
    .avg('duration')
    .avg('build_date')
    .first('flaking')
    .add(DG.STR_AGG.CONCAT_UNIQUE, 'status')
    .add(DG.STR_AGG.CONCAT_UNIQUE, 'result')
    .aggregate();

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

  await pivot.columns.addNewCalculated('stable', generateCommonValueFormula("And", ' concat unique(status)}', '== "passed"'), 'bool');
  await pivot.columns.addNewCalculated('failing', generateCommonValueFormula("Or", ' concat unique(status)}', '== "failed"'), 'bool');
  await pivot.columns.addNewCalculated('flaking', generateCommonValueFormula("Or", ' first(flaking)}', '== true'), 'bool');
  let schemas = await grok.dapi.stickyMeta.getSchemas();
  let schema = schemas.filter((schema) => schema.name == 'Autotests').at(0);
  var meta = await grok.dapi.stickyMeta.getAllValues(schema, pivot.columns.byName('test'));
  pivot.columns.add(meta.col('ignore?'));
  pivot.columns.add(meta.col('ignoreReason'));
  await pivot.columns.addNewCalculated('needs_attention', "And(Or(${flaking},${failing}), Not(${ignore?}))", 'bool');

  function replaceColumn(prefix, type, buildName, newType) {
    var col = pivot.columns.byName(prefix + type);
    col?.setTag('friendlyName', buildName + newType);
  }

  pivot.columns.byName('owner').semType = 'User';
  pivot.columns.byName('owner').setTag('cell.renderer', 'User');
  pivot.columns.byName('test').semType = 'test';

  for (var i = 1; i <= builds.length; i++) {
    var buildName = buildNames[i - 1];
    replaceColumn(i, ' concat unique(status)', i, '');
    replaceColumn(i, ' concat unique(result)', buildName, ' result');
    var colResult = pivot.columns.byName(i + ' concat unique(result)');
    if (colResult !== null)
      colResult.semType = 'stackTrace';
    replaceColumn(i, ' avg(build_date)', buildName, ' build_date');
    replaceColumn(i, ' avg(duration)', buildName, ' duration');
  }


  var attachedTickets = meta.col('tickets');
  attachedTickets.name = 'jira';
  pivot.columns.add(attachedTickets);


  var ticketColumns = 0;
  for (var i = 0; i < pivot.rowCount; i++) {
    var tickets = attachedTickets.get(i)?.split(',') ?? [];
    for (var j = 0; j < tickets.length; j++) {
      pivot.columns.getOrCreate('ticket ' + j, DG.TYPE.STRING).set(i, tickets[j]);
      ticketColumns = Math.max(ticketColumns, j + 1);
    }
  }
  for (var i = 0; i < ticketColumns; i++)
    await pivot.columns.addNewCalculated('severity ' + i, 'JiraConnect:getJiraField(${ticket ' + i + '}, "priority:name")', DG.TYPE.STRING);
  priorityOrders = ['Highest', 'High', 'Medium', 'Low', 'Lowest', '']
  for (var i = 0; i < pivot.rowCount; i++) {
    var maxPriority = 5;
    for (var j = 0; j < ticketColumns; j++) {
      var priority = 0;
      while (priority + 1 < priorityOrders.length && pivot.col('severity ' + j).get(i) != priorityOrders[priority])
        priority++;
      maxPriority = Math.min(priority, maxPriority);
    }
    if (maxPriority >= priorityOrders.length)
      maxPriority = priorityOrders.length - 1;
    pivot.columns.getOrCreate('severity', DG.TYPE.STRING).set(i, priorityOrders[maxPriority]);
  }

  for (var i = 0; i < ticketColumns; i++) {
    pivot.columns.remove('severity ' + i);
    pivot.columns.remove('ticket ' + i);
  }

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

  pivot.name = '0. Tests Dashboard';
  // setTimeout(() => {
  //   grok.functions.call('UsageAnalysis:testDashboard');
  // }, 2000);

  out = pivot;
}



if (result.rowCount == 0)  {
  out = DG.DataFrame.fromColumns([
    DG.Column.fromType('string', 'test', 1),
    DG.Column.fromType('string', 'owner', 1)
  ]);
}
else
  await postprocess();