# MultiPlot

MultiPlot is a [package](https://datagrok.ai/help/develop/develop#packages) for the [Datagrok](https://datagrok.ai)
platform.

MultiPlot was developed in purpose to provide to user ability to easily create dashboards based on complex data which
have common influencing parameter (usually presented as axis X on charts). Common axis X combined with compact form of
presenting charts allow user quickly estimate data meaning.

Distribution of the vertical space among the plots is based on following principles:

1. There are fixed sized plots, size of which specified in pixels, for instance '50px';
2. There are proportional sized plots, size of which speicfied in percents, for instance '25%'
3. There are flex-style space distribution, for instance '2flex'.
4. Each plot could have title. Title could have default size specified in options or individual size.

Types of charts:

1. Scatterplot;
2. Line chart;
3. Bar chart;
4. Status markers (shape and color of markers are used for display data status)

## Examples

To make these examples working please open tables first. Feel free to use tables from 'data' folder or upload your own
tables.

### 1. simple 2 plots

Source code:

```
//name: MultiSimpleTwoPlots
//description: Hello world script
//language: javascript

let options = {
    series: [
        {
            tableName: 'test',
            title: 'plot 0',
            type: 'scatter',
            x: 'col1',
            y: 'col2',
            height: '1flex',
            show: 1,
        },
        {
            tableName: 'test',
            title: 'plot 1',
            type: 'line',
            x: 'col1',
            y: 'col3',
            color: 'red',
            height: '2flex',
            show: 1,
        },
    ]
}

async function func1() {
  let tablesIter = grok.shell.tables;
  let tables = {};
  for (table of tablesIter) {
    tables[table.name] = table;
  };
  let startTable = tables['test'];
  let view = grok.shell.addTableView(startTable);
  let viewer = await startTable.plot.fromType('MultiPlot', {
      paramOptions: JSON.stringify(options),
  });
  view.addViewer(viewer);
}

func1();
```

#### 1. The result of script execution

![image of script result](img/simple_2_plots.png?raw=true "The result")

### 2. advanced example (full features)

This example show following features:

1. Status chart, divided by categories using column and display markers colored by values and limits
2. Timelines chart, display continuous data
3. Multiline chart. Automatically created several linecharts. Categories are taken from specified column and used to
   distribute data through linecharts

<details>
  <summary>Show scriptc</summary>

```
//name: LoadMultiPlot
//description: Hello world script
//language: javascript

let options = {
  series: [
    {
      tableName: 'lb',
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
      tableName: 'ae',
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
      tableName: 'lb',
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
    if (table.name == 'lb') {
      table.filter.init(e => {
        let row = table.row(e);
        return row['USUBJID'] === myId;
      })
    } // lb
  } // tables

  let startTable = tables['lb']
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

```

</details>

#### 2. The result of script execution

![image of script result](img/full_plots_anim.gif?raw=true "The result")
