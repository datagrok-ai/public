//language: javascript
//input: dataframe result
//output: dataframe out

async function postprocess() {
  let batchIndices = result.col('batch_index').categories;
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
      statusCol.colors.setCategorical({
        'passed': '#2ca02c',
        'skipped': '#ffa500',
        'failed': '#8e342a',
      })
    let col = pivot.col(`${batchName} concat unique(result)`);
    for (var j = 0; j < col.length; j++) {
      let tickets = col.getString(j);
      col.set(j, tickets.replaceAll('\n', ','));
    }
    let ticketsStatusCol = await pivot.columns.addNewCalculated(`${batchNames[i]} tickets status`, `UsageAnalysis:getTicketsVerdict(\${${batchNames[i]} concat unique(result)})`, DG.TYPE.STRING);
    if (ticketsStatusCol)
      ticketsStatusCol.colors.setCategorical({
        'Fixed': '#2ca02c',
        'Particialy Fixed': '#ffa500',
        'Wasn\'t Fixed': '#8e342a',
      })
    order.push(`${batchNames[i]} concat unique(result)`);
    order.push(`${batchNames[i]} tickets status`);
  }
 pivot.columns.setOrder(order)

  // let schemas = await grok.dapi.stickyMeta.getSchemas();
  // let schema = schemas.filter((schema) => schema.name == 'Autotests').at(0);

  // var attachedTickets = (await grok.dapi.stickyMeta.getAllValues(schema, pivot.columns.byName('test'))).col('tickets');
  // attachedTickets.name = 'jira';

  // var ticketColumns = 5;
  // for (var i = 0; i < pivot.rowCount; i++) {
  //   var tickets = attachedTickets.get(i)?.split(',') ?? [];
  //   for (var j = 0; j < batchIndices.length; j++)
  //     tickets.push(...pivot.col(`${batchIndices[j]} concat unique(result)`).get(i).split(','));

  //   tickets = [...new Set(tickets)];
  //   var ptr = 0;
  //   for (var j = 0; j < tickets.length; j++) {
  //     var ticketMatch = tickets[j].match(/GROK-\d+/);
  //     if (ticketMatch) {
  //         pivot.columns.getOrCreate('ticket ' + ptr, DG.TYPE.STRING).set(i, ticketMatch[0]);
  //         ptr++;
  //     }
  //     ticketColumns = Math.max(ticketColumns, ptr + 1);
  //   }
  // }

  // for (var i = 0; i < ticketColumns; i++)
  //   await pivot.columns.addNewCalculated(`severity ${i}`, `JiraConnect:getJiraField(RegExpExtract(\${ticket ${i}}, \'GROK-\d+\'), "priority:name")`, DG.TYPE.STRING);

  // priorityOrders = ['Highest', 'High', 'Medium', 'Low', 'Lowest', '']
  // for (var i = 0; i < pivot.rowCount; i++) {
  //   var maxPriority = 5;
  //   for (var j = 0; j < ticketColumns; j++) {
  //     var priority = 0;
  //     while (priority + 1 < priorityOrders.length && pivot.col('severity ' + j).get(i) != priorityOrders[priority])
  //       priority++;
  //     maxPriority = Math.min(priority, maxPriority);
  //   }
  //   if (maxPriority >= priorityOrders.length)
  //     maxPriority = priorityOrders.length - 1;
  //   pivot.columns.getOrCreate('severity', DG.TYPE.STRING).set(i, priorityOrders[maxPriority]);
  // }

  // for (var i = 0; i < ticketColumns; i++) {
  //   pivot.columns.remove('severity ' + i);
  // }

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