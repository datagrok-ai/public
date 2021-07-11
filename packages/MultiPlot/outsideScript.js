//name: Template
//description: Hello world script
//language: javascript

const { color } = require("echarts");

let tables = grok.shell.tables;
for (table of tables) {
  if (table.name == 'ae__2__lb__2_') {
    let view = grok.shell.addTableView(table);
	  let viewer = view.addViewer('MultiPlot', {
          paramA: "string outside", 
          "Param A": "string outside2"
    });
  }
}

// options format proposal:
let options = {
  colors: [

  ],
  titleHeight: '25px',
  timeLineLineWidth: 1,
  timeLineCircleRange: 3,

  plots: [
    {
      table: table1,
      fields: ['field1', 'field2', ...], // default is numeric[0, 1]
      type: 'line', // default is 'scatter'
      axis: [0, 1],   // can be [0, [1,2]] for timelines, default is [0, 1]
      height: '20%', // default is '1flex'
      title: 'title1', // can be omitted, if so then plots become closer to each other
      show: false, // default is true,
      markerColor: (value, i, array) => color(),
      markerType: (value, i, array) => color()
    },
    {
      // ...
    }
  ] // plots
}

// only series:
let series = [{
    table: 'table1',
    type: 'line',
    x: 'date',
    y: 'value',
    color: 'severity',
    marker: 'status',
    height: '20px'
  }];