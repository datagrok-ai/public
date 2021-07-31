//name: LoadMultiPlot
//description: Hello world script
//language: javascript

let options = {
  id: '01-701-1015',
  series: [
    {
      tableName: 'lb2',
      title: 'outside title 0',
      type: 'scatter',
      x: 'LBDY',
      y: 'LBTEST',
      yType: 'category',
      markerShape: 'square',
      height: '1flex',
      show: 1,
      visualMap: {
        type: 'piecewise',
        column: 'LBSEQ',
        pieces: [
          { min: 20, max: 250, color: ['red'] },
        ],
        dimension: 2,
      },
    },
    {
      tableName: 'lb2',
      type: 'scatter',
      x: 'LBDY',
      y: 'LBTEST',
      height: '1flex',
      show: 1,
    },

    // timeLines
    {
      tableName: 'ae__2__lb__2_',
      title: 'outside title 1',
      type: 'timeLine',
      y: 'AETERM',
      x: ['AESTDY', 'AEENDY'],
      yType: 'category',
      color: 'red',
      markerShape: 'square',
      height: '2flex',
      show: 1,
    },

    // linechart with filter 
    {
      tableName: 'lb2',
      title: 'outside title 0',
      type: 'line',
      x: 'LBDY',
      y: 'LBSTRESN',
      splitByColumnName: 'LBTEST',
      maxLimit: 5,
      condition2: {
        field: 'USUBJID',
        value: '01-701-1023'
      },
      yType: 'value',
      markerShape: 'square',
      height: '1flex',
      show: 1,
    },
    /*
            {
        tableName: 'lb2',
        //title: 'outside title 0',
        type: 'line',
        x: 'LBDY',
        y: 'LBSEQ',
        condition2: {
          field: 'USUBJID',
          value: '01-701-1023'
        },
        yType: 'value',
        markerShape: 'square',
        height: '1flex',
        show: 1,
      }
      */
  ]

}

let myId = '01-701-1146'
myId = '01-701-1015'
async function func1() {
  let tablesIter = grok.shell.tables;
  let tables = {}
  for (table of tablesIter) {
    tables[table.name] = table;
  }
  for (table of tablesIter) {
    console.log('Table found: ', table.name);
    if (table.name == 'ae__2__lb__2_') {
      table.filter.init(e => {
        let row = table.row(e);
        return row["USUBJID"] === myId;
      })
    }
    if (table.name == 'ae2') {
      table.filter.init(e => {
        let row = table.row(e);
        return row['USUBJID'] === myId;
      })
    } // ae2
    if (table.name == 'lb2') {
      table.filter.init(e => {
        let row = table.row(e);
        return row['USUBJID'] === myId;
      })
    } // lb2
  } // tables

  let startTable = tables['lb2']
  console.log(startTable)
  let view = grok.shell.addTableView(startTable);
  let viewer = await startTable.plot.fromType('MultiPlot', {
    paramA: "string outside",
    paramOptions: JSON.stringify(options),
  })
  view.addViewer(viewer);
}

func1();
