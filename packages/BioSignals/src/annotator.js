/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

import * as echarts from 'echarts';

export class Annotator extends DG.JsViewer {
    constructor() {
        super();

        let my_data = [[0, 9], [1, 40], [2, 100], [3, 20], [4, 30], [5, 10], [6, 30], [7, 10]];
        let option = { xAxis: {}, yAxis: {}, series: [{data: my_data, type: 'line'}]};

        let chartDiv = ui.div(null, { style: { position: 'absolute', left: '0', right: '0', top: '0', bottom: '0'}} );
        this.root.appendChild(chartDiv);
        this.chart = echarts.init(chartDiv);
        this.chart.setOption(option);
        this.subs.push(ui.onSizeChanged(chartDiv).subscribe((_) => this.chart.resize()));
    }
}