# MultiPlot

MultiPlot is a [package](https://datagrok.ai/help/develop/develop#packages) for the [Datagrok](https://datagrok.ai) platform.

MultiPlot was developed in purpose to provide to user ability to easily create dashboards based on complex data which have common influencing parameter (usually presented as axis X on charts). Common axis X combined with compact form of presenting charts allow user quickly estimate data meaning.

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

##EXAMPLES

To make these examples working please open tables first. Feel free to use tables from 'data' folder or upload your own tables.

1. Simple 2 plots

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

###The result of script execution:

![image of script result](img/simple_2_plots.png?raw=true "The result")
