import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as echarts from 'echarts';


// Based on this example: https://echarts.apache.org/examples/en/editor.html?c=radar
export class RadarViewer extends DG.JsViewer {
  myChart: echarts.ECharts;
  id: number | undefined;

  constructor(id?: number | undefined) {
    super();
    const chartDiv = ui.div([], { style: { position: 'absolute', left: '0', right: '0', top: '0', bottom: '0'}} );
    this.root.appendChild(chartDiv);
    this.myChart = echarts.init(chartDiv);
    this.subs.push(ui.onSizeChanged(chartDiv).subscribe((_) => this.myChart.resize()));
    this.id = id;
  }

  onTableAttached() {
    this.subs.push(this.dataFrame.selection.onChanged.subscribe((_) => this.render()));
    this.subs.push(this.dataFrame.filter.onChanged.subscribe((_) => this.render()));

    this.render();
  }

  render() {
    const option: { [key: string]: any } = {
      radar: {
        name: {
          textStyle: {
            color: '#fff',
            backgroundColor: '#999',
            borderRadius: 3,
            padding: [3, 5],
          },
        },
        indicator: [],
      },
      series: [{
        type: 'radar',
        data: [],
      }],
    };

    const columns = Array.from(this.dataFrame.columns.numerical);
    let data = option.series[0].data;

    if (typeof this.id !== 'undefined') {
      for (const c of columns) {
        option.radar.indicator.push({name: c.name, avg: this.avg(c) });
      }
      this.pushData(data, columns, this.id);
    } else {
      for (const c of columns) {
        option.radar.indicator.push({name: c.name, max: c.max });
      }
      for (let i = 0; i < this.dataFrame.rowCount; i++) {
        this.pushData(data, columns, i);
      }
    }
    this.myChart.setOption(option);
  }

  avg(col: DG.Column) {
    let sum = 0;
    let data = col.getRawData();
    let rowCount = this.dataFrame.rowCount;
    for (let i = 0; i < rowCount; i++) {
      sum += data[i];
    }
  return sum;
  }

  pushData(data: any, columns: DG.Column<any>[], i: number) {
    data.push({
      name: `row ${i}`,
      value: columns.map((c) => c.get(i)),
    });
  }
}

