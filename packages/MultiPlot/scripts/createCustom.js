//name: Template
//description: viewer pass params
//language: javascript

let options = {
    series: [
        {
            table: 'ae__2__lb__2_',
            title: 'outside title 0',
            type: 'scatter',
            x: 'LBDY',
            y: ['LBTEST', 'AEENDY'],
            yType: 'category',
            markerShape: 'square',
            height: '1flex',
            show: 1,
            visualMap: {
                type: 'piecewise',
                pieces: [
                    { min: 20, max: 250, color: ['red'] },
                ],
                dimension: 2,
            },
        },

        // timeLines
        {
            table: 'ae__2__lb__2_',
            title: 'outside title 1',
            type: 'timeLine',
            x: 'LBTEST',
            y: ['AESTDY', 'AEENDY'],
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