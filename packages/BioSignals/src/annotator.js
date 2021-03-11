/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

import * as echarts from 'echarts';

export class Annotator extends DG.JsViewer {
    constructor() {
        super();

        let chartDiv = ui.div(null, { style: { position: 'absolute', left: '0', right: '0', top: '0', bottom: '0'}} );
        this.root.appendChild(chartDiv);
        this.chart = echarts.init(chartDiv);
        this.subs.push(ui.onSizeChanged(chartDiv).subscribe((_) => this.chart.resize()));
    }

    onTableAttached() {
        this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
        this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_) => this.render(false)));

        let columnToPlot = this.dataFrame.columns.byName('testEcg');
        let arrayToPlot = new Array(columnToPlot.length);
        for (let i = 0; i < columnToPlot.length; i++) {
            arrayToPlot[i] = [i, columnToPlot.get(i)];
        }
        let option = {
            xAxis: {},
            yAxis: {},
            series: [
                {
                    data: arrayToPlot,
                    type: 'line',
                    symbolSize: 1,
                    lineStyle: {
                        width: 0.5,
                        type: "dotted"
                    }
                }
            ]
        };

        this.chart.setOption(option);
        this.render();
    }

    detach() {
        this.subs.forEach(sub => sub.unsubscribe());
    }

}