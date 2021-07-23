//name: Template
//description: viewer pass params
//language: javascript

let options = {
    series: [
      {
        table: 'ae__2__lb__2_',
        title: 'outside title1',
        type: 'scatter',
        x: 'AESTDY',
        //    y: 'LBTEST',
        y: 'LBSTRESN',
  
        yType: 'value',
        color: 'red',
        markerShape: 'square',
        height: '1flex',
        show: 1,
     //   itemStyle: {
       //   color: (e) => e.data[2] > 33 ? 'green' : 'defaultColor',
  //      },
        symbol: (e) => e[2] > 22 ? 'rect' : 'circle'
      },
      {
        table: 'ae__2__lb__2_',
        title: 'outside title2',
        type: 'timeLine',
        x: 'LBTEST',
        y: ['AESTDY', 'LBDY'],
        yType: 'category',
        color: 'red',
        markerShape: 'square',
        height: '2flex',
        show: 1,
      },
    ]
  }
  
  async function func1() {
    let tables = grok.shell.tables;
    for (table of tables) {
      if (table.name == 'ae__2__lb__2_') {
        let view = grok.shell.addTableView(table);
        let viewer = await table.plot.fromType('MultiPlot', {
              paramA: "string outside", 
              paramOptions: JSON.stringify(options),
              "Param A": "string outside2"
        })
        view.addViewer(viewer);
      }
    }
    }
    
    func1();