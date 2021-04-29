/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

import * as echarts from 'echarts';

function getTime(now) {
  if (now.getSeconds() < 10)
    return now.getUTCFullYear() + '-' + (now.getUTCMonth() + 1) + '-' + now.getUTCDate() + 'T' + now.getUTCHours() + ':' +
      now.getUTCMinutes() + ':0' + now.getUTCSeconds() + '.' + now.getUTCMilliseconds() + 'Z';
  return now.getUTCFullYear() + '-' + (now.getUTCMonth() + 1) + '-' + now.getUTCDate() + 'T' + now.getUTCHours() + ':' +
    now.getUTCMinutes() + ':' + now.getUTCSeconds() + '.' + now.getUTCMilliseconds() + 'Z';
}

let base = +new Date('Jan 01 1970 00:00:00 GMT+0000');

export class AnnotatorViewer extends DG.JsViewer {
  constructor() {
    super();

    this.chartDiv = ui.div(null, {style: {position: 'absolute', left: '0', right: '0', top: '0', bottom: '0'}});
    this.root.appendChild(this.chartDiv);
    this.chart = echarts.init(this.chartDiv);
    this.subs.push(ui.onSizeChanged(this.chartDiv).subscribe((_) => this.chart.resize()));
  }

  onTableAttached() {
    let signalValues = this.dataFrame.columns.byIndex(0);
    const samplingPeriodInMilliseconds = 1000 / signalValues.getTag('samplingFrequency');
    this.timeOfFirstSample = new Date('Jan 01 1970 00:00:00 GMT+0000');

    let data = new Array(signalValues.length);
    let base1 = base;
    for (let i = 0; i < signalValues.length; i++) {
      let now = new Date(base1 += samplingPeriodInMilliseconds);
      data[i] = [
        getTime(now),
        signalValues.get(i)
      ];
    }
    if (this.dataFrame.columns.byName('indicesOfRPeak')) {
      let indicesOfRPeak = this.dataFrame.columns.byIndex(1);
      let numberOfPointAnnotations = indicesOfRPeak.stats.valueCount;
      indicesOfRPeak = indicesOfRPeak.getRawData();
      indicesOfRPeak = indicesOfRPeak.slice(indicesOfRPeak.length - numberOfPointAnnotations);

      this.markedPoints = new Array(numberOfPointAnnotations);
      let now = new Date(base += indicesOfRPeak[0] * samplingPeriodInMilliseconds);
      for (let i = 1; i < numberOfPointAnnotations; i++) {
        now = new Date(base += (indicesOfRPeak[i] - indicesOfRPeak[i - 1]) * samplingPeriodInMilliseconds);
        this.markedPoints[i] = [
          getTime(now),
          signalValues.get(indicesOfRPeak[i])
        ];
      }
    } else {
      this.markedPoints = [];
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
        text: 'Input (' + signalValues.getTag('samplingFrequency') + ' Hz)',
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
        boundaryGap: false,
        axisLabel: {
          formatter: (function(value) {
            value = new Date(value);
            return (value.getSeconds() < 10) ? value.getMinutes() + ":0" + value.getSeconds() : value.getMinutes() + ":" + value.getSeconds();
          })
        }
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
    this.render(samplingPeriodInMilliseconds);
  }

  render(samplingPeriodInMilliseconds) {
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
        valuesInSearchSegment.push(segment[i][1]);
      }
      let indexOfMax = getIndexOfMaximumValue(valuesInSearchSegment);
      let now = new Date(timeOfFirstSample.getTime() + (indexToSearchFrom + indexOfMax + 1) * samplingPeriod);
      return getTime(now);
    }

    let deltaToSearchMaximum = 10;
    this.chartDiv.addEventListener('click', event => {
      let pointInPixel = [event.offsetX, event.offsetY];
      if (this.chart.containPixel('grid', pointInPixel)) {
        let pointInGrid = this.chart.convertFromPixel('grid', pointInPixel);
        let dateString = getTimeOfMaximumValueNearClick(samplingPeriodInMilliseconds, this.data, pointInGrid[0], deltaToSearchMaximum, this.timeOfFirstSample);
        this.markedPoints.push([dateString, this.data.find((c) => c[0] === dateString)[1]]);
        this.chart.setOption({
          series: [
            {
              data: this.data,
              type: 'line'
            },
            {
              type: 'scatter',
              symbolSize: 10,
              symbol: 'circle',
              data: this.markedPoints
            }
          ]
        })
      }
    })
  }
}