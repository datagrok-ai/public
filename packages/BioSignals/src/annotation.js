import * as echarts from 'echarts';

let chartDom = document.getElementById('main');
let myChart = echarts.init(chartDom);
let option;

function indexOfMin(arr) {
    if (arr.length === 0) {return -1;}
    let min = arr[0];
    let minIndex = 0;
    for (let i = 0; i < arr.length; i++) {
        if (arr[i] < min) {
            minIndex = i;
            min = arr[i];
        }
    }
    return minIndex;
}

function getClosestPointIndex(point, arr) {
    if (arr.length === 0) {
        return 0;
    }
    let a = [];
    for (let i = 0; i < arr.length; i++) {
        let v = Math.abs( point - arr[i][0] );
        a.push( v );
    }
    return indexOfMin(a);
}

let my_data = [[0,9], [1, 40], [2, 100], [3, 20], [4, 30], [5, 10], [6, 30], [7, 10]];
let markedPoints = [];

option = {
    xAxis: {},
    yAxis: {},
    series: [
        {
            data: my_data,
            type: 'line'
        }
    ]
};
let zr = myChart.getZr();
myChart.setOption(option);
zr.on('dblclick', function (params) {
    let pointInPixel = [params.offsetX, params.offsetY];
    let pointInGrid = myChart.convertFromPixel('grid', pointInPixel);
    if (myChart.containPixel('grid', pointInPixel)) {
        let ind = Math.round(pointInGrid[0]);
        markedPoints.push([ind, my_data[ind][1]]);
        myChart.setOption({
            yAxis: {},
            series: [
                {
                    data: my_data,
                    type: 'line'
                },
                {
                    type: 'scatter',
                    symbolSize: 20,
                    symbol: 'circle',
                    data: markedPoints
                }
            ]
        })
    }
});

zr.on('contextmenu', function (params) {
    var pointInPixel = [params.offsetX, params.offsetY];
    var pointInGrid = myChart.convertFromPixel('grid', pointInPixel);
    if (myChart.containPixel('grid', pointInPixel)) {
        let ind = getClosestPointIndex(pointInGrid[0], markedPoints);
        console.log('ind1 = ' + ind);
        console.log('markedPoints.length1 = ' + markedPoints.length);
        markedPoints.splice(ind, 1);
        console.log('markedPoints.length2 = ' + markedPoints.length);
        myChart.setOption({
            yAxis: {},
            series: [
                {
                    data: my_data,
                    type: 'line'
                },
                {
                    type: 'scatter',
                    symbolSize: 20,
                    symbol: 'circle',
                    data: markedPoints
                }
            ]
        })
    }
});

option && myChart.setOption(option);