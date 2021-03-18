/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

import * as echarts from 'echarts';

export class AnnotatorViewer extends DG.JsViewer {
  constructor() {
    super();

    this.chartDiv = ui.div(null, { style: { position: 'absolute', left: '0', right: '0', top: '0', bottom: '0'}} );
    this.root.appendChild(this.chartDiv);
    this.chart = echarts.init(this.chartDiv);
    this.subs.push(ui.onSizeChanged(this.chartDiv).subscribe((_) => this.chart.resize()));
  }

  onTableAttached() {
    let signalValues = this.dataFrame.columns.byName('testEcg');
    let base = +new Date();
    this.timeOfFirstSample = new Date();
    this.samplingPeriod = 24 * 3600 * 1000;

    let data = new Array(signalValues.length);
    for (let i = 0; i < signalValues.length; i++) {
      let now = new Date(base += this.samplingPeriod);
      data[i] = [
        [now.getFullYear(), now.getMonth() + 1, now.getDate()].join('/'),
        signalValues.get(i)
      ];
    }

    let indicesOfRPeak = this.dataFrame.columns.byName('indicesOfRPeak');
    let numberOfPointAnnotations = indicesOfRPeak.stats.valueCount;
    indicesOfRPeak = indicesOfRPeak.getRawData();
    indicesOfRPeak = indicesOfRPeak.slice(indicesOfRPeak.length - numberOfPointAnnotations);

    this.markedPoints = new Array(numberOfPointAnnotations);
    base = +new Date();
    let samplingPeriod = 24 * 3600 * 1000;
    let now = new Date(base += indicesOfRPeak[0] * samplingPeriod);
    for (let i = 1; i < numberOfPointAnnotations; i++) {
      now = new Date(base += (indicesOfRPeak[i] - indicesOfRPeak[i-1]) * samplingPeriod);
      this.markedPoints[i] = [
        [now.getFullYear(), now.getMonth() + 1, now.getDate()].join('/'),
        signalValues.get(indicesOfRPeak[i])
      ];
    }

    this.data = data;

    this.option = {
      tooltip: {
        trigger: 'axis',
        position: function (pt) {
          return [pt[0], '10%'];
        }
      },
      title: {
        left: 'center',
        text: 'Input signal',
      },
      toolbox: {
        feature: {
          dataZoom: {
            yAxisIndex: 'none'
          },
          restore: {},
          saveAsImage: {}
        }
      },
      xAxis: {
        type: 'time',
        boundaryGap: false
      },
      yAxis: {
        type: 'value',
        boundaryGap: [0, '100%']
      },
      dataZoom: [
        {
          type: 'inside',
          start: 0,
          end: 10
        },
        {
          start: 0,
          end: 10
        }
      ],
      series: [
        {
          type: 'line',
          symbol: 'none',
          data: data
        },
        {
          type: 'scatter',
          symbolSize: 20,
          symbol: 'circle',
          data: this.markedPoints
        }
      ]
    };

    this.chart.setOption(this.option);
    this.render();
  }

  render() {
    function getDateStringFromUnixTimestamp(unix_timestamp) {
      let now = new Date(unix_timestamp);
      return [now.getFullYear(), now.getMonth() + 1, now.getDate()].join('/');
    }

    function getIndexOfMaximumValue(arr) {
      let max = arr[0];
      let maxIndex = 0;
      for (let i = 0; i < arr.length; i++) {
        if (arr[i] > max) {
          maxIndex = i;
          max = arr[i];
        }
      }
      return maxIndex;
    }

    function getSampleIndex(unix_timestamp, timeOfFirstSample, samplingPeriod) {
      return Math.round((unix_timestamp - timeOfFirstSample) / samplingPeriod);
    }

    function getTimeOfMaximumValueNearClick(samplingPeriod, data, unixTimestampOfClick, searchDelta, timeOfFirstSample) {
      let clickedIndex = getSampleIndex(unixTimestampOfClick, timeOfFirstSample, samplingPeriod);
      let indexToSearchFrom = clickedIndex - searchDelta;
      let indexToSearchTo = clickedIndex + searchDelta;
      let segment = data.slice(indexToSearchFrom, indexToSearchTo);
      let valuesInSearchSegment = [];
      for (let i = 0; i < segment.length; i++) {
        valuesInSearchSegment.push( segment[i][1] );
      }
      let indexOfMax = getIndexOfMaximumValue(valuesInSearchSegment);
      let now = new Date(timeOfFirstSample.getTime() + (indexToSearchFrom + indexOfMax + 1) * samplingPeriod);
      return [now.getFullYear(), now.getMonth() + 1, now.getDate()].join('/');
    }

    let deltaToSearchMaximum = 10;
    this.chartDiv.addEventListener('click', event => {
      let pointInPixel = [event.offsetX, event.offsetY];
      if (this.chart.containPixel('grid', pointInPixel)) {
        let pointInGrid = this.chart.convertFromPixel('grid', pointInPixel);
        let dateString = getTimeOfMaximumValueNearClick(this.samplingPeriod, this.data, pointInGrid[0], deltaToSearchMaximum, this.timeOfFirstSample);
        this.markedPoints.push([dateString, this.data.find((c) => c[0] === dateString)[1]]);
        this.chart.setOption({
          series: [
            {
              data: this.data,
              type: 'line'
            },
            {
              type: 'scatter',
              symbolSize: 20,
              symbol: 'circle',
              data: this.markedPoints
            }
          ]
        })
      }
    })
  }
}