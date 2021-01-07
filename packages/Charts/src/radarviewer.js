import * as echarts from 'echarts';


// Based on this example: https://echarts.apache.org/examples/en/editor.html?c=radar
export class RadarViewer extends DG.JsViewer {
    constructor() {
        super();
        let chartDiv = ui.div(null, { style: { position: 'absolute', left: '0', right: '0', top: '0', bottom: '0'}} );
        this.root.appendChild(chartDiv);
        this.myChart = echarts.init(chartDiv);
        this.subs.push(ui.onSizeChanged(chartDiv).subscribe((_) => this.myChart.resize()));
    }

    onTableAttached() {
        this.subs.push(this.dataFrame.selection.onChanged.subscribe((_) => this.render()));
        this.subs.push(this.dataFrame.filter.onChanged.subscribe((_) => this.render()));

        this.render();
    }

    render() {

        let option = {
            radar: {
                name: {
                    textStyle: {
                        color: '#fff',
                        backgroundColor: '#999',
                        borderRadius: 3,
                        padding: [3, 5]
                    }
                },
                indicator: [ ]
            },
            series: [{
                type: 'radar',
                data: [ ]
            }]
        };

        let columns = Array.from(this.dataFrame.columns.numerical);

        for (let c of columns)
          option.radar.indicator.push( {name: c.name, max: c.max });

        let data = option.series[0].data;
        for (let i = 0; i < this.dataFrame.rowCount; i++) {
            data.push({
                name: `row ${i}`,
                value: columns.map((c) => c.get(i))
            });
        }

        this.myChart.setOption(option);
    }
}
