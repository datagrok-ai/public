//language: javascript
//input: dataframe result
//output: dataframe out

async function postprocess() {
  let builds = result.col('build_index').categories;
  let buildNames = result.col('build').categories.sort();

  let pivot = result
    .groupBy(['test'])
    .pivot('build_index')
    // .avg('build_date')
    .add(DG.STR_AGG.CONCAT_UNIQUE, 'status')
    // .add(DG.STR_AGG.CONCAT_UNIQUE, 'result')
    .aggregate();



  function replaceColumn(prefix, type, buildName, newType) {
    var col = pivot.columns.byName(prefix + type);
    col?.setTag('friendlyName', buildName + newType);
    if (col !== null) 
      col.name = buildName;
  }

  pivot.columns.byName('test').semType = 'autotest';

  for (var i = 0; i < builds.length; i++) {
    var buildName = buildNames[i];
    replaceColumn(builds[i], ' concat unique(status)', (i + 1).toString(), '');
    // replaceColumn(i, ' concat unique(result)', buildName, ' result');
    // replaceColumn(i, ' avg(build_date)', buildName, ' build_date');
  }

  let schemas = await grok.dapi.stickyMeta.getSchemas();
  let schema = schemas.filter((schema) => schema.name == 'Autotests').at(0);

  var attachedTickets = (await grok.dapi.stickyMeta.getAllValues(schema, pivot.columns.byName('test'))).col('tickets');
  attachedTickets.name = 'jira';
  // pivot.columns.add(attachedTickets);


  var ticketColumns = 5;
  for (var i = 0; i < pivot.rowCount; i++) {
    var tickets = attachedTickets.get(i)?.split(',') ?? [];
    for (var j = 0; j < tickets.length; j++)   {
      pivot.columns.getOrCreate('ticket ' + j, DG.TYPE.STRING).set(i, tickets[j]);
      ticketColumns = Math.max(ticketColumns, j + 1);
    }
  }
  for (var i = 0; i < ticketColumns; i++)
    await pivot.columns.addNewCalculated(`severity ${i}`, `JiraConnect:getJiraField(RegExpExtract(\${ticket ${i}}, \'GROK-\d+\'), "priority:name")`, DG.TYPE.STRING);
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
    // pivot.columns.remove('ticket ' + i);
  }
  // console.log(pivot);


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
  pivot.name = 'test track dashboard';

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