//name: LoadMultiPlot
//description: Hello world script
//language: javascript

let options = {
  series: [
    {
      tableName: 'lb2',
      title: 'statusChart',
      type: 'scatter',
      x: 'LBDY',
      y: 'LBTEST',
      // extraFields is an array to load into echart data arrays
      // all fields later combined into one array [x, y, extraFields]
      // user can address fields by index, for instance index 3 means field "LBORNRLO"
      extraFields: ['LBORRES', 'LBORNRLO', 'LBORNRHI'],
      yType: 'category',                // can be 'value' or 'category'
      statusChart: {
        valueField: 2,                  // index of field with test value
        splitByColumnName: 'LBTEST',    // column to get categories
        categories: ["Basophils", "Urate", "Glucose"], // fixed categories
        minField: 3,                    // min and max normal value 
        maxField: 4,                    // will be displayed with default color, otherwises "red"
        maxLimit: 5,                    // max number of categories
        alertColor: 'red',
      },
      markerShape: 'circle',
      height: '1flex',                  // height can be '30px', '20%', '3flex'
      show: 1,
    },

    // timeLines
    {
      tableName: 'ae__2__lb__2_',
      title: 'Timelines',
      type: 'timeLine',
      y: 'AETERM',                      // category column
      x: ['AESTDY', 'AEENDY'],          // [startTime, endTime]
      yType: 'category',
      color: 'red',                     // color of marker
      markerShape: 'circle',
      height: '2flex',
      show: 1,
    },

    // multi linechart 
    {
      tableName: 'lb2',
      title: 'Multi Linechat',
      type: 'line',
      multi: true,
      x: 'LBDY',
      y: 'LBSTRESN',
      splitByColumnName: 'LBTEST',                    // get categories from this column
      categories: ["Basophils", "Urate", "Glucose"],  // fixed categories
      maxLimit: 5,                                    // max number of linecharts 
      yType: 'value',
      markerShape: 'square',
      height: '1flex',
      show: 1,
    },
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
        return 1; // no filter to show all lines
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
    paramOptions: JSON.stringify(options),
  })
  setTimeout((e) => {
    viewer.setOptions({
      testField1: 'testValue1',
    });
  }, 5000)
  view.addViewer(viewer);
}

func1();
