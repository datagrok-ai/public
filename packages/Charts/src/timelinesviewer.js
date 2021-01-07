import * as echarts from 'echarts';

export class TimelinesViewer extends DG.JsViewer {
    constructor() {
        super();
        
//        setTimeout(function() { 

            console.log('init');
            let chartDiv = ui.div(null, { style: { position: 'absolute', left: '0', right: '0', top: '0', bottom: '0'}} );
            this.root.appendChild(chartDiv);
            this.myChart = echarts.init(chartDiv);
            ui.onSizeChanged(chartDiv).subscribe((_) => this.myChart.resize());
    
            var option = {
                title: {
                    text: 'ECharts entry example'
                },
                tooltip: {},
                legend: {
                    data:['Sales']
                },
                xAxis: {
                    data: ["shirt","cardign","chiffon shirt","pants","heels","socks"]
                },
                yAxis: {},
                series: [{
                    name: 'Sales',
                    type: 'bar',
                    data: [5, 20, 36, 10, 10, 20]
                }]
            };
    
            // use configuration item and data specified to show chart
            this.myChart.setOption(option);

//        }.bind(this), 50);
    }

    // init() {
    //     this.render();
    // }

    // onTableAttached() {
    //     // console.log('attached');
    //     // console.log(this.dataFrame);
    //     // debugger;
    //     // this.subs.push(this.dataFrame.selection.onChanged.subscribe((_) => this.render()));
    //     // this.subs.push(this.dataFrame.filter.onChanged.subscribe((_) => this.render()));

    //     this.render();
    // }

}
